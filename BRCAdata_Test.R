################################################################################

# load packages
library(r.jive)
library(tidyverse)
library(GenomicDataCommons)
library(missMDA)
library(dplyr)
library(impute)

# load data
data(BRCA_data)
mRNA <- Data$miRNA
methyl <- Data$Methylation

# Set some columns to missing at random (in one data set e.g. mRNA)
set.seed(500) # For reproducibility
missing_fraction <- 0.1
mask <- matrix(runif(nrow(mRNA) * ncol(mRNA)) < missing_fraction, nrow(mRNA), ncol(mRNA))
mRNA_missing <- mRNA
mRNA_missing[mask] <- NA


##############################################INITIALIZE VALUES FOR EM ALGORITHM
impute.init <- function(x) {
  # 1) Initialize using row and column averages
  mindex <- is.na(x)
  mcol <- colMeans(x, na.rm=T)
  mrow <- rowMeans(x, na.rm=T)
  mtot <- mean(x, na.rm=T)
  last <- 0
  xtemp <- x
  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)) {
      # For single missing values, average col and row means
      xtemp[i,j] <- (mrow[i] + mcol[j]) / 2
    }
  }
  for (i in 1:nrow(x)) {
    if (is.na(mrow[i])) {
      xtemp[i,] <- mcol
    }
  }
  for (j in 1:ncol(x)) {
    if (is.na(mcol[j])) {
      xtemp[,j] <- mrow
    }
  }
  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)) {
      if (is.na(mcol[j]) & is.na(mrow[i])) {
        xtemp[i,j] <- mtot
      }
    }
  }
  xtemp[!mindex] <- x[!mindex]
  
  return(list(xtemp = xtemp, mindex = mindex)) #return(xtemp)
}

# PCA
impute.pca <- function (x, maxiter=500) {
  
  #impute starting values, get transpose of data matrix
  # xtemp <- impute.init(as.matrix(x))
  # dimnames(xtemp) <- NULL
  init_result <- impute.init(as.matrix(x))
  xtemp <- init_result$xtemp
  mindex <- init_result$mindex
  dimnames(xtemp) <- NULL
  
  conv <- F
  niter <- 0
  convlist <- c()
  #this is the part that needs changed I think 
  while(!conv & niter < maxiter) {
    
    # 2) get loadings and scors with PCA, estimate new x with this equation
    pca <- prcomp(xtemp, scale=F, center=F) # HERE IS THE CHANGE
    loadings <- pca$rotation
    scores <- pca$x
    xfull <- loadings %*% t(scores)
    
    # 3) Impute missing values using PCA in (2) (EM algorithm)
    last <- xtemp
    xtemp[mindex] <- xfull[mindex]
    
    # Check for convergence (leaving this the same for now)
    # get the norm of the matrix - "f" specifies the Frobenius norm
    #     (the Euclidean norm of x treated as if it were a vector);
    convlist[niter+1] <- norm(xtemp - last, type='f')^2 / norm(xtemp, type='f')^2
    if (norm(xtemp - last, type='f')^2 / norm(xtemp, type='f')^2 < .0001) { conv <- T }
    print(paste("Iteration ", niter,": ",convlist[niter+1],sep=""))
    niter <- niter + 1
  }
  
  return(list(x=xtemp,nconv=niter,convlist=convlist))
}



# Use the impute.pca function to impute the missing values in mRNA_missing
result_pca <- impute.pca(mRNA_missing)

# Extract the imputed dataset
mRNA_imputed_pca <- result_pca$x

#################################################################################################

# R.JIVE

impute.init <- function(x) {
  x <- as.matrix(x)
  x[is.na(x)] <- 0
  return(x)
}

#x and y need to start as matrices
impute.jive <- function(x, y, maxiter=500) {
  
  # sdx <- sd(as.matrix(x), na.rm=T)
  # sdy <- sd(as.matrix(y), na.rm=T)
  
  #impute starting values, get transpose of data matrices
  xtemp <- impute.init(x)
  dimnames(xtemp) <- NULL
  mx_index <- is.na(x)
  #xtemp <- xtemp/sdx #NEW
  
  ytemp <- impute.init(y)
  dimnames(ytemp) <- NULL
  my_index <- is.na(y)
  # ytemp <- impute.init(as.matrix(y))
  # dimnames(ytemp) <- NULL
  # my_index <- is.na(y)
  #ytemp <- ytemp/sdy #NEW
  
  conv <- F
  niter <- 0
  convlistx <- c()
  convlisty <- c()
  
  while(!conv & niter < maxiter) {
    
    Data <- list(xtemp, ytemp)
    
    # 2) run JIVE - if the equations below are correct
    results = jive(Data, rankJ=1, rankA=c(1,3), method="given",scale = T, center = T)
    #results = jive(Data, scale = T, center = T)
    xfull = results$joint[[1]] + results$indiv[[1]]
    yfull = results$joint[[2]] + results$indiv[[2]]
    
    xfull = xfull*results$scale$`Scale Values`[1]
    yfull = yfull*results$scale$`Scale Values`[2]
    
    for (r in 1:nrow(xfull)) {
      xfull[r,] = xfull[r,] + results$scale$`Center Values`[[1]][r]
    }
    for (r in 1:nrow(yfull)) {
      yfull[r,] = yfull[r,] + results$scale$`Center Values`[[2]][r]
    }
    
    # 3) Impute missing values using PCA in (2) (EM algorithm)
    lastx <- xtemp
    lasty <- ytemp
    xtemp[mx_index] <- xfull[mx_index]
    ytemp[my_index] <- yfull[my_index]
    
    # Check for convergence (leaving this the same for now)
    # get the norm of the matrix - "f" specifies the Frobenius norm
    #     (the Euclidean norm of x treated as if it were a vector);
    convlistx[niter+1] <- norm(xtemp - lastx, type='f')^2 / norm(xtemp, type='f')^2
    convlisty[niter+1] <- norm(ytemp - lasty, type='f')^2 / norm(ytemp, type='f')^2
    
    if ((norm(xtemp - lastx, type='f')^2 / norm(xtemp, type='f')^2 < .0001) &
        (norm(ytemp - lasty, type='f')^2 / norm(ytemp, type='f')^2 < .0001)) { conv <- T }
    print(paste("Iteration ", niter,": ",convlistx[niter+1]," ,",convlisty[niter+1], sep=""))
    niter <- niter + 1
    
  }
  return(list(x=xtemp,y=ytemp,nconv=niter,convlistx=convlistx, convlisty=convlisty))
}

# Use the impute.jive function to impute the missing values in mRNA_missing
result_jive <- impute.jive(x = as.data.frame(mRNA_missing), y = as.data.frame(methyl))


# Extract the imputed dataset
mRNA_imputed_jive <- result_jive$x


# Output csv file
write.csv(mRNA_missing, file = "mRNA_missing.csv")
write.csv(mRNA_imputed_pca, file = "mRNA_imputed_pca.csv")
write.csv(mRNA_imputed_jive, file = "mRNA_imputed_jive.csv")

# Read results
setwd("...")
mRNA_pca <- read.csv("mRNA_imputed_pca.csv")[,-1]
mRNA_jive <- read.csv("mRNA_imputed_jive.csv")[,-1]

# RMSE
sqrt(mean((as.matrix(mRNA_pca) - mRNA)^2))
sqrt(mean((as.matrix(mRNA_jive) - mRNA)^2))

# Correlation
cor(as.vector(as.matrix(mRNA_pca)), as.vector(mRNA))
cor(as.vector(as.matrix(mRNA_jive)), as.vector(mRNA))

# Row_means
mRNA_means <- rowMeans(mRNA) %*% t(rep(1,348))
mRNA_rmi <- mRNA
mRNA_rmi[mask] <- mRNA_means[mask]

# RMSE
sqrt(mean((mRNA_rmi - mRNA)^2))




##########################################################################################################
##### Test the full column with NA result ###################################################################
##########################################################################################################

set.seed(500) # For reproducibility
mRNA_misscol <- mRNA
mRNA_misscol[,"TCGA.A1.A0SJ.01A.11R"] <- NA

##############################################INITIALIZE VALUES FOR EM ALGORITHM
impute.init <- function(x) {
  # 1) Initialize using row and column averages
  mindex <- is.na(x)
  mcol <- colMeans(x, na.rm=T)
  mrow <- rowMeans(x, na.rm=T)
  mtot <- mean(x, na.rm=T)
  last <- 0
  xtemp <- x
  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)) {
      # For single missing values, average col and row means
      xtemp[i,j] <- (mrow[i] + mcol[j]) / 2
    }
  }
  for (i in 1:nrow(x)) {
    if (is.na(mrow[i])) {
      xtemp[i,] <- mcol
    }
  }
  for (j in 1:ncol(x)) {
    if (is.na(mcol[j])) {
      xtemp[,j] <- mrow
    }
  }
  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)) {
      if (is.na(mcol[j]) & is.na(mrow[i])) {
        xtemp[i,j] <- mtot
      }
    }
  }
  xtemp[!mindex] <- x[!mindex]
  
  return(list(xtemp = xtemp, mindex = mindex)) #return(xtemp)
}

# PCA
impute.pca <- function (x, maxiter=500) {
  
  #impute starting values, get transpose of data matrix
  # xtemp <- impute.init(as.matrix(x))
  # dimnames(xtemp) <- NULL
  init_result <- impute.init(as.matrix(x))
  xtemp <- init_result$xtemp
  mindex <- init_result$mindex
  dimnames(xtemp) <- NULL
  
  conv <- F
  niter <- 0
  convlist <- c()
  #this is the part that needs changed I think 
  while(!conv & niter < maxiter) {
    
    # 2) get loadings and scors with PCA, estimate new x with this equation
    pca <- prcomp(xtemp, scale=F, center=F) # HERE IS THE CHANGE
    loadings <- pca$rotation
    scores <- pca$x
    xfull <- loadings %*% t(scores)
    
    # 3) Impute missing values using PCA in (2) (EM algorithm)
    last <- xtemp
    xtemp[mindex] <- xfull[mindex]
    
    # Check for convergence (leaving this the same for now)
    # get the norm of the matrix - "f" specifies the Frobenius norm
    #     (the Euclidean norm of x treated as if it were a vector);
    convlist[niter+1] <- norm(xtemp - last, type='f')^2 / norm(xtemp, type='f')^2
    if (norm(xtemp - last, type='f')^2 / norm(xtemp, type='f')^2 < .0001) { conv <- T }
    print(paste("Iteration ", niter,": ",convlist[niter+1],sep=""))
    niter <- niter + 1
  }
  
  return(list(x=xtemp,nconv=niter,convlist=convlist))
}

# Use the new impute.pca function to impute the missing values in mRNA_missing
col_result_pca <- impute.pca(mRNA_misscol)
# Extract the imputed dataset
mRNA_col_imputed_pca <- col_result_pca$x


# Use the impute.jive function to impute the missing values in mRNA_missing
col_result_jive <- impute.jive(x = as.data.frame(mRNA_misscol), y = as.data.frame(methyl))
# Extract the imputed dataset
mRNA_col_imputed_jive <- row_result_jive$x
