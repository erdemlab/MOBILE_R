##### MOBILE -omics integrator (Code by Cemal Erdem, Ph.D.)
##### Last update: June 2022
##### Integration of transcriptomic and epigenomic datasets 
#from matlab to R

setwd(".../Mobile_R")

# #Data/matrix pre-processing

RNGkind(kind = "Mersenne-Twister")
set.seed(6)
XXraw <- read.table("RNAseq_lvl42_data.txt", header = FALSE, colClasses = rep("numeric", 15))
YYraw <- read.table("RPPA_lvl4_data.txt", header = FALSE, sep = ",",colClasses = rep("numeric", 15))
std_n <- function(x) {
  x <- as.numeric(x)  # Ensure x is numeric
  sqrt(sum((x - mean(x, na.rm = TRUE))^2, na.rm = TRUE) / length(x))
}
std_XXraw <- apply(XXraw, 1, std_n)
std_YYraw <- apply(YYraw, 1, std_n)
XX_indices <- 1:nrow(XXraw)
YY_indices <- 1:nrow(YYraw)
XX_combined <- cbind(std_XXraw, XX_indices)
YY_combined <- cbind(std_YYraw, YY_indices)
XX_sorted <- XX_combined[order(-XX_combined[, 1]), ]
YY_sorted <- YY_combined[order(-YY_combined[, 1]), ]
XXvar <- XX_sorted
YYvar <- YY_sorted
XXindices2keep <- floor(0.1 * nrow(XXraw))
YYindices2keep <- floor(0.2 * nrow(YYraw))
XX_indices <- order(XXvar[1:XXindices2keep, 2])
YY_indices <- order(YYvar[1:YYindices2keep, 2])
XX_sorted_indices <- XXvar[1:XXindices2keep, 2][XX_indices]
YY_sorted_indices <- YYvar[1:YYindices2keep, 2][YY_indices]
XX <- XXraw[XX_sorted_indices, , drop = FALSE]
YY <- YYraw[YY_sorted_indices, , drop = FALSE]
RNAseqIDs <- as.matrix(XXvar[1:XXindices2keep, 2][order(XXvar[1:XXindices2keep, 2])])
RPPAIDs <- as.matrix(YYvar[1:YYindices2keep, 2][order(YYvar[1:YYindices2keep, 2])])

# define the prepnormats function
# Function to process/normalize columns and rows of matrices
# By Cemal Erdem, Ph.D. (June 2022)
prepnormmats <- function(X, rc, B) {
  # Convert to matrix if X is a data frame and remove row/column names
  X <- as.matrix(X)
  rownames(X) <- NULL
  colnames(X) <- NULL
  
  rowsnum <- nrow(X)
  colmnsnum <- ncol(X)
  Xn <- X
  
  if (rc == 1) {  # row center
    for (qq in 1:rowsnum) {
      muj <- mean(Xn[qq, ], na.rm = TRUE)
      Xn[qq, ] <- Xn[qq, ] - muj
    }
  } else if (rc == 2) {  # row len = 1
    for (qq in 1:rowsnum) {
      norm2j <- sqrt(sum(Xn[qq, ]^2, na.rm = TRUE))
      Xn[qq, ] <- Xn[qq, ] / norm2j
    }
  } else if (rc == 3) {  # column center
    for (qq in 1:colmnsnum) {
      muj <- mean(Xn[, qq], na.rm = TRUE)
      Xn[, qq] <- Xn[, qq] - muj
    }
  } else if (rc == 4) {  # column len = 1
    for (qq in 1:colmnsnum) {
      norm2j <- sqrt(sum(Xn[, qq]^2, na.rm = TRUE))
      Xn[, qq] <- Xn[, qq] / norm2j
    }
  } else if (rc == 5) {  # mat len = A
    A0 <- sqrt(sum(X^2, na.rm = TRUE))
    Xn <- sqrt(B) * (X / A0)
  } else if (rc == 6) {  # column-center, row-center, row-len=1
    for (qq in 1:colmnsnum) {
      muj <- mean(Xn[, qq], na.rm = TRUE)
      Xn[, qq] <- Xn[, qq] - muj
    }
    
    for (qq in 1:rowsnum) {
      muj <- mean(Xn[qq, ], na.rm = TRUE)
      Xn[qq, ] <- Xn[qq, ] - muj
      norm2j <- sqrt(sum(Xn[qq, ]^2, na.rm = TRUE))
      Xn[qq, ] <- Xn[qq, ] / norm2j
    }
  } else if (rc == 7) {  # column-center, row-center, row-std=1
    for (qq in 1:colmnsnum) {
      muj <- mean(Xn[, qq], na.rm = TRUE)
      Xn[, qq] <- Xn[, qq] - muj
    }
    
    for (qq in 1:rowsnum) {
      muj <- mean(Xn[qq, ], na.rm = TRUE)
      Xn[qq, ] <- Xn[qq, ] - muj
      stdj <- sd(Xn[qq, ], na.rm = TRUE)
      Xn[qq, ] <- Xn[qq, ] / stdj
    }
  }
  
  return(Xn)
}

XX3 <- prepnormmats(XX, 6, 1)
YY3 <- prepnormmats(YY, 6, 1)

# Load necessary libraries
library(glmnet)
library(parallel)
library(caret)

# Set parameters
toplot <- 0
rtimes <- 10  # Number of Lasso modules to run. Use 2 to verify installation, set to 10000 for MOBILE
CVnfolds <- 4

# Create a new directory for outputs
dir.create("FULL_RpR", showWarnings = FALSE)

# define the glmnetRunner function
glmnetRunner <- function(XX, YY, toplot = 0, CVnfolds = 3) {
  library(glmnet)
  XX <- as.matrix(XX)
  YY <- as.matrix(YY)
  if (ncol(XX) != ncol(YY)) {
    stop("Number of columns in XX must match the number of columns in YY.")
  }
  nosamples <- nrow(YY)
  nfeatures <- nrow(XX)
  TTg <- matrix(0, nrow = nosamples, ncol = nfeatures)
  a0Vals <- numeric(nosamples)
  count2 <- 0
  for (qq in 1:nosamples) {
    y_vector <- YY[qq, ]
    BBg <- cv.glmnet(t(XX), y_vector, family = "gaussian", type.measure = "mse", nfolds = CVnfolds, parallel = TRUE, standardize = FALSE)
    Lming <- BBg$lambda.min
    Lindexg <- which(BBg$lambda == Lming)
    Bmin <- BBg$glmnet.fit$beta[, Lindexg]
    if (any(Bmin != 0)) {
      TTg[qq, ] <- as.vector(Bmin)
      a0Vals[qq] <- BBg$glmnet.fit$a0[Lindexg]
    } else {
      Lming2 <- BBg$lambda.1se
      Lindexg2 <- which(BBg$lambda == Lming2)
      Bmin2 <- BBg$glmnet.fit$beta[, Lindexg2]
      if (any(Bmin2 != 0)) {
        TTg[qq, ] <- as.vector(Bmin2)
        a0Vals[qq] <- BBg$glmnet.fit$a0[Lindexg2]
      } else {
        count2 <- count2 + 1
        DFuniqs <- unique(BBg$glmnet.fit$df)
        DFreps <- table(BBg$glmnet.fit$df)
        Dfmaxid <- which.max(DFreps)
        Dfmaxid2 <- DFuniqs[Dfmaxid]
        Dfmaxid3 <- which(BBg$glmnet.fit$df == Dfmaxid2)
        Bmaxid <- Dfmaxid3[length(Dfmaxid3)]
        Bmin3 <- BBg$glmnet.fit$beta[, Bmaxid]
        TTg[qq, ] <- as.vector(Bmin3)
        a0Vals[qq] <- BBg$glmnet.fit$a0[Bmaxid]
      }
    }
  }
  devg2 <- TTg %*% XX - YY
  dsedevg2 <- sum(devg2^2)
  YYpg <- TTg %*% XX
  
  # Plotting if required
  if (toplot == 1) {
    plot(1:nrow(YY), YY[, 1], type = "l", col = "blue", lwd = 2, ylim = range(c(YY[, 1], YYpg[, 1])), xlab = "Samples", ylab = "Values")
    lines(1:nrow(YY), YYpg[, 1], col = "red", lwd = 2, lty = 2)
  }
  return(list(TTg = TTg, BBg = BBg, a0Vals = a0Vals, YYpg = YYpg, dsedevg2 = dsedevg2))
}



# Define the parsavef function 
parsavef <- function(fname, TLas1, BLas1, a0Vals1, YYpg, YYpgdev) {
  save(TLas1, BLas1, a0Vals1, YYpg, YYpgdev, file = fname)
}

# Define the parclearf function 
parclearf <- function() {
  list(TLas1 = NULL, BLas1 = NULL, a0Vals1 = NULL, YYpg = NULL, YYpgdev = NULL)
}


# Run the simulations
cl <- makeCluster(detectCores() - 1)
clusterEvalQ(cl, {
  library(glmnet)
  library(caret)
})

clusterExport(cl, c("glmnetRunner", "parsavef", "parclearf", "XX3", "YY3", "toplot", "CVnfolds"))

parLapply(cl, 1:rtimes, function(rr) {
  rndgntr <- sample.int(1e8, 1)
  set.seed(rndgntr)
  result <- glmnetRunner(XX3, YY3, toplot, CVnfolds)
  fname <- paste0("FULL_RpR/TLas_r", rr, ".RData")
  parsavef(fname, result$TTg, result$BBg, result$a0Vals, result$YYpg, result$dsedevg2)
  parclearf()
})

stopCluster(cl)

## Simulation clean-up and save Lasso Coefficient Matrices

# Load necessary libraries
library(Matrix)

# Set directory and get list of files
dir_path <- "FULL_RpR"
listTLasFiles <- list.files(path = dir_path, pattern = "TLas_r.*\\.RData$", full.names = TRUE)
numFiles <- length(listTLasFiles)
TLasMats <- vector("list", numFiles)
qqtodel <- vector()

# Loop for saving coefficient matrices
for (qq in 1:numFiles) {
  file_info <- file.info(listTLasFiles[qq])
  if (file_info$size > 1000) {
    name1 <- listTLasFiles[qq]
    load(name1)  # This loads the RData file, which should contain TLas1 and/or TLas2
    if (exists("TLas1")) {
      if (!is.matrix(TLas1)) {
        TLas1 <- as.matrix(TLas1)
      }
      TLas1_dense <- as.matrix(TLas1)
      NumName1 <- as.numeric(gsub(".*_r(\\d+)\\.RData$", "\\1", name1))
      TLasMats[[qq]] <- list(matrix = TLas1_dense, NumName1 = NumName1)
      
      rm(TLas1)
    }
    if (exists("TLas2")) {
      if (!is.matrix(TLas2)) {
        TLas2 <- as.matrix(TLas2)
      }
      TLas2_dense <- as.matrix(TLas2)
      NumName2 <- as.numeric(gsub(".*_r(\\d+)\\.RData$", "\\1", name1))
      TLasMats[[qq + numFiles]] <- list(matrix = TLas2_dense, NumName2 = NumName2)
      rm(TLas2)
    }
    cat(qq, "\n")
  }
}

# Create a list to store y-intercept value arrays
a0Mats <- vector("list", numFiles)

# Loop for saving y-intercept value arrays (usually all zero)
for (qq in 1:numFiles) {
  file_info <- file.info(listTLasFiles[qq])
  if (file_info$size > 1000) {
    name1 <- listTLasFiles[qq]
    load(name1)  # This loads the RData file, which should contain a0Vals1
    if (exists("a0Vals1")) {
      if (!is.matrix(a0Vals1)) {
        a0Vals1 <- as.matrix(a0Vals1)
      }
      aname1 <- gregexpr("_", name1)
      NumName1 <- as.numeric(substr(name1, aname1[[1]][length(aname1[[1]])] + 2, nchar(name1) - 4))
      a0Mats[[qq]] <- list(matrix = a0Vals1, NumName1 = NumName1)
      rm(a0Vals1)
    }
    cat(qq, "\n")
  }
}

# Remove elements in qqtodel from the lists
if (length(qqtodel) > 0) {
  listTLasFiles <- listTLasFiles[-qqtodel]
  a0Mats <- a0Mats[-qqtodel]
  TLasMats <- TLasMats[-qqtodel]
}

# Alternatively, save all data in a single .rds file
save_data <- list(listTLasFiles = listTLasFiles, a0Mats = a0Mats, TLasMats = TLasMats, qqtodel = qqtodel)
saveRDS(save_data, file = "FULL_RpR/RpR_FULL.rds")

## Finding Robust Lasso Coefficient Matrix

runAssocRankerRpR <- function(cutoffval, TLasMats, XX_Data, XX_IDs, YY_Data, YY_IDs) {
  matrices <- lapply(TLasMats, function(x) x$matrix)
  numnetworks <- length(matrices)
  numedges <- numeric(numnetworks)
  TLasAllsum_1s <- matrix(0, nrow = nrow(matrices[[1]]), ncol = ncol(matrices[[1]]))
  TLasAllsum <- matrix(0, nrow = nrow(matrices[[1]]), ncol = ncol(matrices[[1]]))
  for (qq in seq_len(numnetworks)) {
    TLastemp <- matrices[[qq]]
    numedges[qq] <- sum(TLastemp != 0) 
    TLasAllsum <- TLasAllsum + TLastemp
    TLastemp[TLastemp != 0] <- 1 
    TLasAllsum_1s <- TLasAllsum_1s + TLastemp
  }
  cutofflimit <- floor(numnetworks * cutoffval)
  TLas <- (TLasAllsum_1s >= cutofflimit)
  TLasAllmean <- TLasAllsum / TLasAllsum_1s
  TLasAllmean[is.na(TLasAllmean)] <- 0
  TLasTarget <- TLas * TLasAllmean
  TLasTarget[is.na(TLasTarget)] <- 0
  idx <- which(TLasTarget != 0, arr.ind = TRUE)
  intranks <- data.frame(
    jj = idx[, 2],
    ii = idx[, 1],
    ensembl_id_XX = as.character(XX_Data$ensembl_id[XX_IDs[idx[, 2]]]),
    hgnc_id_XX = as.character(XX_Data$hgnc_id[XX_IDs[idx[, 2]]]),
    empty1 = rep(" ", length(idx[, 1])),
    empty2 = rep(" ", length(idx[, 1])),
    hgnc_id_YY = as.character(YY_Data$hgnc_id[YY_IDs[idx[, 1]]]),
    empty3 = rep(" ", length(idx[, 1])),
    score = TLasTarget[idx],
    abs_score = abs(TLasTarget[idx]),
    antibodyfull = as.character(YY_Data$antibodyfull[YY_IDs[idx[, 1]]]),
    pathway = as.character(YY_Data$pathway[YY_IDs[idx[, 1]]])
  )
  intranks_r <- intranks[order(-intranks$abs_score), ]
  is_na_or_empty <- is.na(intranks_r$hgnc_id_XX) | intranks_r$hgnc_id_XX == "" | intranks_r$hgnc_id_XX == "NA"
  intranks_r$hgnc_id_XX[is_na_or_empty] <- intranks_r$ensembl_id_XX[is_na_or_empty]
  intranks_ii <- intranks_r[intranks_r$hgnc_id_XX == intranks_r$hgnc_id_YY, ]
 return(list(intranks_r, intranks_ii, numedges, TLasTarget))
}
RPPA<-readRDS(".../Mobile_R/RPPA.rds")
RNAseq_lvl42<-readRDS(".../Mobile_R/RNAseq_lvl42.rds")
TLasMats<-readRDS(".../Mobile_R/FULL_RpR/TLasMats.rds")
TLasMats<-as.matrix(TLasMats)
results <- runAssocRankerRpR(0.5, TLasMats, RNAseq_lvl42, RNAseqIDs,RPPA,RPPAIDs)
intranksR_RpR_FULL <- results[[1]]
intsii_RpR_FULL <- results[[2]]
numedges_RpR_FULL <- as.matrix(results[[3]])
TLasTarget_RpR_FULL <- results[[4]]
# Change the below name according to inputs: For example, for EGF-LOGO -> TLasBestRA_EGFLOGO
saveRDS(list(intranksR_RpR_FULL = intranksR_RpR_FULL, intsii_RpR_FULL = intsii_RpR_FULL, numedges_RpR_FULL = numedges_RpR_FULL, TLasTarget_RpR_FULL = TLasTarget_RpR_FULL), file = "FULL_RpR/TLasBestRpR_FULL.rds")
cat("DONE and saved\n")


