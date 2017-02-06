# Random fourier transformation
randomFourierTrans <- function (X, Dim, sigma, seedW=NULL) {
    # Draw weights of multivariate normal distribution
    n <- dim(X) [1]
    d <- dim(X) [2]
    X <- t(X)
    set.seed(seedW)
    rW <- rmvnorm (n=d, sigma=diag(Dim) / sigma)
    
    # Calculate Z matrix
    inner <- crossprod(rW, X)
    Z <- rbind(cos(inner), sin(inner)) / sqrt(Dim)
    Output <- list(Z=Z, rW=rW)
    return(Output)
}

# Used for prediction of new inputs given weight matrix
fourierTransPredict <- function (newx, rW) {
  newx <- t(newx)
  inner <- crossprod(rW, newx)
  Z <- rbind(cos(inner), sin(inner)) / sqrt(dim(rW)[2])
  return(Z)
}

# Robust Standardization
# X: Design Matrix
# Standardizes a matrix with median and median absolute deviation
robustStandard <- function (X) {
  colIndex <- 1:dim(X) [2]
  MedianVal <- sapply(colIndex, function (j) median(x=X [, j]))
  MadValue <- sapply(colIndex, function (j) mad (x=X [, j], constant=1))
  TypeScaling <- ifelse(MadValue!=0, "Robust", "Standard")
  for(j in 1:length(TypeScaling)) {
    if(TypeScaling[j]=="Standard") {
      MedianVal[j] <- mean(X[, j])
      MadValue[j] <- sd(X[, j])
      
      # Only center and set standard deviation to 1
      if(MadValue[j]==0) {
        TypeScaling[j] <- "OnlyCenter"
        MedianVal[j] <- mean(X[, j])
        MadValue[j] <- 1
      }
    }
    X [, j] <- (X [, j] - MedianVal [j]) / MadValue [j]
  }
  attr(X, which="StandValues") <- list(Location=MedianVal, 
                                       Scale=MadValue, 
                                       Type=TypeScaling)
  return(X)
}

# Fit KDSN
fitKDSN <- function (y, X, levels=1, 
                     Dim=rep(round(sqrt(dim(X)[1]) / 2), levels), 
                     sigma=rep(1, levels), lambdaRel=rep(0.5, levels), 
                     alpha=rep(0, levels), 
                     info=FALSE, seedW=NULL, standX=TRUE, standY=FALSE,
                     varSelect=rep(FALSE, levels), varRanking=NULL, 
                     varExFreq=rep(NA, levels), 
                     dropHidden=rep(FALSE, levels), 
                     dropHiddenProb=rep(NA, levels), 
                     seedDrop=NULL,
                     baggingInd=NULL,
                     randSubsetInd=NULL,
                     maxItOpt=10000) {
  
  # With variable selection or random subsets
  # Necessary to reset number of columns with dimBaseOrg
  CheckAnyVarSelect <- any(varSelect)
  CheckAnyRandSubsets <- any(sapply(1:length(randSubsetInd), 
                                    function(i) !is.null(randSubsetInd[[i]]) ) )
  if(CheckAnyVarSelect | CheckAnyRandSubsets) {
    if(CheckAnyVarSelect) {
      # Check if varRanking is not NULL
      stopifnot(!is.null(varRanking))
    }
    # Extract number of variables of original X
    dimBaseOrg <- ncol(X)
    dimBase <- dimBaseOrg
    varRankingTemp <- varRanking
  } 

  # Conversion to match format of response and covariates
  Y <- matrix(y, ncol=1)
  
  # Checks of Inputs
  stopifnot (is.matrix(Y) & is.matrix(X))
  stopifnot (dim(Y)[2]==1)
  stopifnot (length(sigma) == levels)
  stopifnot (length(lambdaRel) == levels)
  stopifnot (length(alpha) == levels)
  stopifnot (length(seedW) == levels || is.null(seedW))
  
  # Preparation
  rftX <- vector("list", levels)
  linWeights <- vector("list", levels)
  StandValues <- list(Design=list(Location=rep(NA, dim(X)[2]), 
                                  Scale=rep(NA, dim(X)[2]), 
                                  Type=rep(NA, dim(X)[2])), 
                      LinPreds=list(Location=rep(NA, levels-1), 
                                    Scale=rep(NA, levels-1), 
                                    Type=rep(NA, levels-1)))
  StandValuesY <- list(Location=NA, Scale=NA, Type=NA)
  
  # Standardization with median and median absolute deviation
  if(standX) {
    X <- robustStandard (X=X)
    StandValues$Design <- attr(X, which="StandValues")
  }
  
  # Standardization of responses
  if(standY) {
    Y <- robustStandard (X=Y)
    StandValuesY <- attr(Y, which="StandValues")
  }

  ################
  # Only one level
  
  if(levels==1) {
    
    # Bagging
    if(!is.null(baggingInd[[1]])) {
      X <- X[baggingInd[[1]], , drop=FALSE]
    }
    
    # Random subsets
    if(!is.null(randSubsetInd[[1]])) {
      X <- X[, randSubsetInd[[1]] , drop=FALSE]
      if(varSelect[1]) {
        dimBase <- length(randSubsetInd[[1]])
        varRankingTemp <- varRanking[varRanking %in% randSubsetInd[[1]] ]
      }
    }

    # Variable selection
    if(varSelect[1]) {
      
      # Only original covariates may be removed
      if(round(dimBase*varExFreq[1]) >= 1) {
        
        # Exclude variables
        exCol <- varRankingTemp[1:round(dimBase*varExFreq[1])]
        
        # If it is the first level and all variables should be dropped,
        # then include only the best variable, 
        # because otherwise the algorithm is instable
        if(length(exCol) == length(varRankingTemp)) {
          if(ncol(X)==1) {
            exCol <- 2
          }
          else{
            exCol <- exCol[-length(exCol)]
          }
        }
        X_new <- X [, -exCol, drop=FALSE]
        
        # Apply random fourier transformation
        rftX [[1]] <- randomFourierTrans (X=X_new, Dim=Dim[1], 
                                          sigma=sigma[1], seedW=seedW[1])
      }
      else{
        # Apply random fourier transformation
        rftX [[1]] <- randomFourierTrans (X=X, Dim=Dim[1], 
                                          sigma=sigma[1], seedW=seedW[1])
      }
    }
    else{
      # Apply random Fourier transformation
      rftX [[1]] <- randomFourierTrans (X=X, Dim=Dim[1], sigma=sigma[1], seedW=seedW[1])
    }
    
    # Calculate y_new and X_new
    y_new <- rftX [[1]] $Z %*% Y
    X_new <- tcrossprod(rftX [[1]] $Z)
    
    # Dropout
    if(dropHidden[1]) {
      set.seed(seedDrop[1])
      randBin <- sample(x=c(1/dropHiddenProb[1], 0), size=prod(dim(X_new)), 
                        prob=c(dropHiddenProb[1], 1-dropHiddenProb[1]),
                        replace=TRUE)
      # To avoid error of zero variance in all predictors in function glmnet
      # use normal distribution noise
      if(sum(randBin)==0) {
        set.seed(seedDrop[1])
        randBin <- rnorm(n=length(randBin), mean=1, 
                         sd=sqrt(dropHiddenProb[1]*(1-dropHiddenProb[1])))
      }
      tempMat <- matrix(randBin, nrow=nrow(X_new), ncol=ncol(X_new))
      X_new <- X_new * tempMat
    }
    
    # Estimate linear weights
    glmnetFit <- glmnet (x=X_new, y=y_new, alpha=alpha[1], 
                         standardize = FALSE, nlambda=1000, maxit=maxItOpt)
    lambdaActual <- lambdaRel [1]*(max(glmnetFit$lambda) - min(glmnetFit$lambda)) + 
      min(glmnetFit$lambda)
    linWeights [[1]] <- coef(glmnetFit, s=lambdaActual, exact=FALSE)
    linWeights [[1]] <- c(as.matrix(linWeights [[1]]))
    
    # Save input
    Input <- list(levels=levels, Dim=Dim, sigma=sigma, 
                  lambdaRel=lambdaRel, alpha=alpha, 
                  info=info, seedW=seedW, standX=standX, standY=standY,
                  varSelect=varSelect, varRanking=varRanking, varExFreq=varExFreq,
                  dropHidden=dropHidden, dropHiddenProb=dropHiddenProb, 
                  seedDrop=seedDrop, baggingInd=baggingInd, randSubsetInd=randSubsetInd)
    
    # Arrange output
    Output <- list(rftX=rftX, linWeights=linWeights, StandValues=StandValues, 
                   StandValuesY=StandValuesY, lambdaActual=lambdaActual)
    All <- list(Output=Output, Input=Input)
    class(All) <- "KDSN"
    return(All)
  }
  
  ###################
  # Two or more levels
  else {
    
    lambdaActual <- vector("numeric", levels)
    for(l in 1:(levels-1)) {
      
      # Bagging
      if(!is.null(baggingInd[[l]])) {
        X_new <- X[baggingInd[[l]], , drop=FALSE]
      }
      
      # Random subsets
      if(!is.null(randSubsetInd[[l]])) {
        
        # Project choosen inclusion indices to equivalent exclusion indices
        if(length(setdiff(1:dimBaseOrg, randSubsetInd[[l]]))!=0){
          newExInd <- -setdiff(1:dimBaseOrg, randSubsetInd[[l]])
        } 
        else{
          newExInd <- -(dimBaseOrg+levels)
        }
        
        if(!is.null(baggingInd[[l]])) {
          X_new <- X_new[, newExInd, drop=FALSE]
        } 
        else{
          X_new <- X[, newExInd, drop=FALSE]
        }
        if(varSelect[l]) {
          dimBase <- length(randSubsetInd[[l]])
          varRankingTemp <- varRanking[varRanking %in% randSubsetInd[[l]] ]
        }
      } 
      else{ # Reset to default, if no random subsets are applied in this level
        if(varSelect[l]) {
          dimBase <- dimBaseOrg
          varRankingTemp <- varRanking
        }
      }
      
      # Variable selection
      if(varSelect[l]) {
        
        # Only original covariates may be removed
        if(round(dimBase*varExFreq[l]) >= 1) {
          
          # Exclude variables
          exCol <- varRankingTemp[1:round(dimBase*varExFreq[l])]
          
          # If it is the first level and all variables should be dropped,
          # then include only the best (last) variable, 
          # because otherwise the algorithm is instable
          if((length(exCol) == length(varRankingTemp)) & l==1) {
            checkLog <- !is.null(baggingInd[[l]]) | !is.null(randSubsetInd[[l]])
            if(checkLog) {
              if(ncol(X_new)==1) {
                exCol <- 2
              }
              else{
                exCol <- exCol[-length(exCol)]
              }
            }
            if(!checkLog) {
              if(ncol(X)==1) {
                exCol <- 2
              }
              else{
                exCol <- exCol[-length(exCol)]
              }
            }
          }
          
          # Take X_new, if either bagging or random subsets are enabled in this level
          if(!is.null(baggingInd[[l]]) | !is.null(randSubsetInd[[l]])) {
            X_new <- X_new [, -exCol, drop=FALSE]
          }
          else{
            X_new <- X [, -exCol, drop=FALSE]
          }
          
          # Apply random fourier transformation
          rftX [[l]] <- randomFourierTrans (X=X_new, Dim=Dim[l], 
                                            sigma=sigma[l], seedW=seedW[l])
        }
        # No variable is selected for exclusion
        else{
          # Take X_new, if either bagging or random subsets are enabled in this level
          if(!is.null(baggingInd[[l]]) | !is.null(randSubsetInd[[l]])) {
            # Apply random fourier transformation
            rftX [[l]] <- randomFourierTrans (X=X_new, Dim=Dim[l], 
                                              sigma=sigma[l], seedW=seedW[l])
          }
          else{
            rftX [[l]] <- randomFourierTrans (X=X, Dim=Dim[l], 
                                              sigma=sigma[l], seedW=seedW[l])
          }
        }
      }
      # Without variable selection
      else{
        if(!is.null(baggingInd[[l]]) | !is.null(randSubsetInd[[l]])) {
          rftX [[l]] <- randomFourierTrans (X=X_new, Dim=Dim [l], 
                                            sigma=sigma [l], seedW=seedW [l])
        }
        else{
          rftX [[l]] <- randomFourierTrans (X=X, Dim=Dim[l], 
                                            sigma=sigma[l], seedW=seedW[l])
        }
      }
      
      # Calculate y_new and X_new
      y_new <- rftX [[l]] $Z %*% Y
      X_new <- tcrossprod(rftX [[l]] $Z)
      
      # Dropout
      if(dropHidden[l]) {
        set.seed(seedDrop[l])
        randBin <- sample(x=c(1/dropHiddenProb[l], 0), size=prod(dim(X_new)), 
                          prob=c(dropHiddenProb[l], 1-dropHiddenProb[l]),
                          replace=TRUE)
        # To avoid error of zero variance in all predictors in function glmnet
        # use normal distribution noise
        if(sum(randBin)==0) {
          set.seed(seedDrop[l])
          randBin <- rnorm(n=length(randBin), mean=1, 
                           sd=sqrt(dropHiddenProb[l]*(1-dropHiddenProb[l])))
        }
        tempMat <- matrix(randBin, nrow=nrow(X_new), ncol=ncol(X_new))
        X_new <- X_new * tempMat
      }
      
      # Estimate linear weights
      glmnetFit <- glmnet (x=X_new, y=y_new, alpha=alpha[l], 
                           standardize = FALSE, nlambda=1000, maxit=maxItOpt)
      lambdaActual[l] <- lambdaRel [l]*(max(glmnetFit$lambda) - min(glmnetFit$lambda)) + 
        min(glmnetFit$lambda)
      linWeights [[l]] <- coef(glmnetFit, s=lambdaActual[l], exact=FALSE)
      linWeights [[l]] <- c(as.matrix(linWeights [[l]])) 
      
      # Predict linear outputs
      linPreds <- cbind(1, t(rftX [[l]] $Z)) %*% linWeights [[l]]
      
      # Expand Input Matrix
      if(standX) {
        StandValues$LinPreds$Location[l] <- median(x=linPreds)
        StandValues$LinPreds$Scale[l] <- mad(x=linPreds, constant=1)
        StandValues$LinPreds$Type[l] <- ifelse(StandValues$LinPreds$Scale[l]!=0, "Robust", "Standard")
        
        if(StandValues$LinPreds$Type[l]=="Standard") {
          StandValues$LinPreds$Location[l] <- mean(linPreds)
          StandValues$LinPreds$Scale[l] <- sd(linPreds)
          
          # Only center if variance is zero
          if(StandValues$LinPreds$Scale[l]==0) {
            StandValues$LinPreds$Type[l] <- "OnlyCenter"
            StandValues$LinPreds$Location[l] <- mean(linPreds)
            StandValues$LinPreds$Scale[l] <- 1
          }
        }
        linPreds <- (linPreds - StandValues$LinPreds$Location[l]) / StandValues$LinPreds$Scale[l]
      }
      X <- cbind(X, linPreds)
      if(info) {cat("level", l, "fit", "\n")}
    }
    
    ################
    # Final level
    
    # Bagging
    if(!is.null(baggingInd[[levels]])) {
      X <- X[baggingInd[[levels]], , drop=FALSE]
    }
    
    # Random subsets
    if(!is.null(randSubsetInd[[levels]])) {
      
      # Project choosen inclusion indices to equivalent exclusion indices
      if(length(setdiff(1:dimBaseOrg, randSubsetInd[[levels]]))!=0){
        newExInd <- -setdiff(1:dimBaseOrg, randSubsetInd[[levels]])
      } 
      else{
        newExInd <- -(dimBaseOrg+levels)
      }
      
      X <- X[, newExInd , drop=FALSE]
      
      # Adapt variable selection dimension and ranking to random subset
      if(varSelect[levels]) {
        dimBase <- length(randSubsetInd[[levels]])
        varRankingTemp <- varRanking[varRanking %in% randSubsetInd[[levels]] ]
      }
    }
    # Reset to default, if no random subsets are applied in this level
    else{
      if(varSelect[levels]) {
        dimBase <- dimBaseOrg
        varRankingTemp <- varRanking
      }
    }
    
    # Variable selection
    if(varSelect[levels]) {
      
      # Only original covariates may be removed
      if(round(dimBase*varExFreq[levels]) >= 1) {
        
        # Exclude variables
        exCol <- varRankingTemp[1:round(dimBase*varExFreq[levels])]
        X_new <- X [, -exCol, drop=FALSE]
        
        # Apply random fourier transformation
        rftX [[levels]] <- randomFourierTrans (X=X_new, Dim=Dim[levels], 
                                          sigma=sigma[levels], seedW=seedW[levels])
      }
      else{
        # Apply random fourier transformation
        rftX [[levels]] <- randomFourierTrans (X=X, Dim=Dim[levels], 
                                          sigma=sigma[levels], seedW=seedW[levels])
      }
    }
    
    # Without variable selection
    else{
      # Apply random Fourier transformation
      rftX [[levels]] <- randomFourierTrans (X=X, Dim=Dim [levels], 
                                             sigma=sigma [levels], seedW=seedW [levels])
    }
    
    # Calculate y_new and X_new
    y_new <- rftX [[levels]] $Z %*% Y
    X_new <- tcrossprod(rftX [[levels]] $Z)
    
    # Dropout
    if(dropHidden[levels]) {
      set.seed(seedDrop[levels])
      randBin <- sample(x=c(1/dropHiddenProb[levels], 0), size=prod(dim(X_new)), 
                        prob=c(dropHiddenProb[levels], 1-dropHiddenProb[levels]),
                        replace=TRUE)
      # To avoid error of zero variance in all predictors in function glmnet
      # use normal distribution noise
      if(sum(randBin)==0) {
        set.seed(seedDrop[levels])
        randBin <- rnorm(n=length(randBin), mean=1, 
                         sd=sqrt(dropHiddenProb[levels]*(1-dropHiddenProb[levels])))
      }
      tempMat <- matrix(randBin, nrow=nrow(X_new), ncol=ncol(X_new))
      X_new <- X_new * tempMat
    }
    
    # Use glmnet as final level
    glmnetFit <- glmnet (x=X_new, y=y_new, alpha=alpha [levels], 
                         standardize = FALSE, nlambda=1000, maxit=maxItOpt)
    lambdaActual[levels] <- lambdaRel [levels]*(max(glmnetFit$lambda) - min(glmnetFit$lambda)) + 
      min(glmnetFit$lambda)
    linWeights [[levels]] <- coef(glmnetFit, s=lambdaActual[levels], exact=FALSE)
    linWeights [[levels]] <- c(as.matrix(linWeights [[levels]]))
    if(info) {cat("level", levels, "fit", "\n")}
    
    # Output
    # Save input
    Input <- list(levels=levels, Dim=Dim, sigma=sigma, lambdaRel=lambdaRel, 
                  alpha=alpha, info=info, seedW=seedW, 
                  standX=standX, standY=standY, varSelect=varSelect, 
                  varRanking=varRanking, varExFreq=varExFreq, dropHidden=dropHidden, 
                  dropHiddenProb=dropHiddenProb, seedDrop=seedDrop,
                  baggingInd=baggingInd, randSubsetInd=randSubsetInd)
    
    # Arrange output
    Output <- list(rftX=rftX, linWeights=linWeights, StandValues=StandValues,
                   StandValuesY=StandValuesY, lambdaActual=lambdaActual)
    All <- list(Output=Output, Input=Input)
    class(All) <- "KDSN"
    return (All)
  }
}

# Predict method for KDSN
predict.KDSN <- function (object, newx, ...) {

    # Variable pre selection
    if(any(names(attributes(object))=="preSelectVars")) {
      newx <- newx[, attr(object, "preSelectVars"), drop=FALSE]
    }
    
    # Variable selection
    if(any(object$Input$varSelect)) {
      # Extract dimension of covariables
      dimBaseOrg <- ncol(newx)
      dimBase <- dimBaseOrg
      varRankingTemp <- object$Input$varRanking
    }
      
    # Preprocessing new data
    noLevels <- object$Input$levels
    if(object$Input$standX) {
      for(j in 1:dim(newx)[2]) {
          newx[, j] <- (newx[, j] - object$Output$StandValues$Design$Location[j]) / 
            object$Output$StandValues$Design$Scale[j]
      }
    }
  
    ############
    # One level
    
    if(noLevels == 1) {
      
      # Random subsets
      if(!is.null(object$Input$randSubsetInd[[1]])) {
        newx <- newx[, object$Input$randSubsetInd[[1]], drop=FALSE]
        if(object$Input$varSelect[1]) {
          dimBase <- length(object$Input$randSubsetInd[[1]])
          varRankingTemp <- object$Input$varRanking[object$Input$varRanking %in% 
                                                      object$Input$randSubsetInd[[1]] ]
        }
      }
      
      # Variable selection
      if(object$Input$varSelect[1]) {
        
        # Only original covariates may be removed
        if(round(dimBase*object$Input$varExFreq[1]) >= 1) {
          
          # Exclude variables
          exCol <- varRankingTemp[1:round(dimBase*object$Input$varExFreq[1])]
          
          # If it is the first level and all variables should be dropped,
          # then include only the best variable, 
          # because otherwise the algorithm is instable
          if(length(exCol) == length(varRankingTemp)) {
            if(ncol(newx)==1) {
              exCol <- 2
            }
            else{
              exCol <- exCol[-length(exCol)]
            }
          }
          newxTemp <- newx [, -exCol, drop=FALSE]
          
          # Transform new data to Fourier space
          tZ <- fourierTransPredict (newx=newxTemp, rW=object$Output$rftX [[1]] $rW)
        }
        else{
          # Transform new data to Fourier space
          tZ <- fourierTransPredict (newx=newx, rW=object$Output$rftX [[1]] $rW)
        }
      }
      else{
        # Transform new data to Fourier space
        tZ <- fourierTransPredict (newx=newx, rW=object$Output$rftX [[1]] $rW)
      }
      
      # Predict linear outputs
      preds <- cbind(1, t(tZ)) %*% object$Output$linWeights [[1]]
      
      # Transform back to original scale, if response has been transformed
      if(object$Input$standY) {
        preds <- preds*object$Output$StandValuesY$Scale + 
          object$Output$StandValuesY$Location
      }
      
      return(c(preds))
    }
    
    #############
    # More levels
  
    else {
      
      # Predict hidden levels
      
      for(l in 1:(noLevels-1)) {
        
        # Random subsets
        if(!is.null(object$Input$randSubsetInd[[l]])) {
          
          # Project choosen inclusion indices to equivalent exclusion indices
          if(length(setdiff(1:dimBaseOrg, object$Input$randSubsetInd[[l]]))!=0){
            newExInd <- -setdiff(1:dimBaseOrg, object$Input$randSubsetInd[[l]])
          } 
          else{
            newExInd <- -(dimBaseOrg + noLevels)
          }
          
          newxTemp <- newx[, newExInd, drop=FALSE]
          if(object$Input$varSelect[l]) {
            dimBase <- length(object$Input$randSubsetInd[[l]])
            varRankingTemp <- object$Input$varRanking[object$Input$varRanking %in% 
                                                        object$Input$randSubsetInd[[l]] ]
          }
        }
        else{
          if(object$Input$varSelect[l]) {
            dimBase <- dimBaseOrg
            varRankingTemp <- object$Input$varRanking
          }
        }
        
        # Variable selection
        if(object$Input$varSelect[l]) {
          
          # Only original covariates may be removed
          if(round(dimBase*object$Input$varExFreq[l]) >= 1) {
            
            # Exclude variables
            exCol <- varRankingTemp[1:round(dimBase*object$Input$varExFreq[l])]
            
            # If it is the first level and all variables should be dropped,
            # then include only the best (last) variable, 
            # because otherwise the algorithm would be instable
            if((length(exCol) == length(varRankingTemp)) & l==1) {
              checkLog <- !is.null(object$Input$randSubsetInd[[l]])
              if(checkLog) {
                if(ncol(newxTemp)==1) {
                  exCol <- 2
                }
                else{
                  exCol <- exCol[-length(exCol)]
                }
              }
              if(!checkLog) {
                if(ncol(newx)==1) {
                  exCol <- 2
                }
                else{
                  exCol <- exCol[-length(exCol)]
                }
              }
            }
            if(!is.null(object$Input$randSubsetInd[[l]])) {
              newxTemp <- newxTemp [, -exCol, drop=FALSE]
            }
            else{
              newxTemp <- newx [, -exCol, drop=FALSE]
            }

            # Transform new data to Fourier space
            tZ <- fourierTransPredict (newx=newxTemp, rW=object$Output$rftX [[l]] $rW)
          }
          else{
            if(!is.null(object$Input$randSubsetInd[[l]])) {
              # Transform new data to Fourier space
              tZ <- fourierTransPredict (newx=newxTemp, rW=object$Output$rftX [[l]] $rW)
            }
            else{
              # Transform new data to Fourier space
              tZ <- fourierTransPredict (newx=newx, rW=object$Output$rftX [[l]] $rW)
            }
          }
        }
        else{
          if(!is.null(object$Input$randSubsetInd[[l]])) {
            # Transform new data to Fourier space
            tZ <- fourierTransPredict (newx=newxTemp, rW=object$Output$rftX [[l]] $rW)
          }
          else{
            # Transform new data to Fourier space
            tZ <- fourierTransPredict (newx=newx, rW=object$Output$rftX [[l]] $rW)
          }
        }
        
        # Predict output values
        preds <- cbind(1, t(tZ)) %*% object$Output$linWeights [[l]]
  
        # Expand inputs
        if(object$Input$standX) {
            preds <- (preds - object$Output$StandValues$LinPreds$Location[l]) / 
              object$Output$StandValues$LinPreds$Scale[l]
        }
        newx <- cbind(newx, preds)
      }
      # Predict output level
      
      # Random subsets
      if(!is.null(object$Input$randSubsetInd[[noLevels]])) {
        
        # Project choosen inclusion indices to equivalent exclusion indices
        if(length(setdiff(1:dimBaseOrg, object$Input$randSubsetInd[[noLevels]]))!=0){
          newExInd <- -setdiff(1:dimBaseOrg, object$Input$randSubsetInd[[noLevels]])
        } 
        else{
          newExInd <- -(dimBaseOrg + noLevels)
        }
        
        newx <- newx[, newExInd, drop=FALSE]
        if(object$Input$varSelect[noLevels]) {
          dimBase <- length(object$Input$randSubsetInd[[noLevels]])
          varRankingTemp <- object$Input$varRanking[object$Input$varRanking %in% 
                                                      object$Input$randSubsetInd[[noLevels]] ]
        }
      }
      else{
        if(object$Input$varSelect[noLevels]) {
          dimBase <- dimBaseOrg
          varRankingTemp <- object$Input$varRanking
        }
      }
      
      # Variable selection
      if(object$Input$varSelect[noLevels]) {
        
        # Only original covariates may be removed
        if(round(dimBase*object$Input$varExFreq[noLevels]) >= 1) {
          
          # Exclude variables
          exCol <- varRankingTemp[1:round(dimBase*object$Input$varExFreq[noLevels])]
          newxTemp <- newx [, -exCol, drop=FALSE]
          
          # Transform new data to Fourier space
          tZ <- fourierTransPredict (newx=newxTemp, rW=object$Output$rftX [[noLevels]] $rW)
        }
        else{
          # Transform new data to Fourier space
          tZ <- fourierTransPredict (newx=newx, rW=object$Output$rftX [[noLevels]] $rW)
        }
      }
      else{
        # Transform new data to Fourier space
        tZ <- fourierTransPredict (newx=newx, rW=object$Output$rftX [[noLevels]] $rW)
      }
      
      # Predict output values
      preds <- cbind(1, t(tZ)) %*% object$Output$linWeights [[noLevels]]
      
      # Transform back to original scale, if response has been transformed
      if(object$Input$standY) {
        preds <- preds*object$Output$StandValuesY$Scale + 
          object$Output$StandValuesY$Location
      }
      
      return(c(preds))
    }
}

######################
# Ensemble predictions

predict.KDSNensemble <- function (object, newx, ...) {
  ensembleSize <- length(object)
  predMat <- matrix(NA, nrow=nrow(newx), ncol=ensembleSize)
  for(i in 1:ensembleSize) {
    varPreSelectInput <- any(names(attributes(object[[i]]))=="preSelectVars")
    if(varPreSelectInput) {
      predMat[, i] <- predict(object=object[[i]], newx=newx)
#                              newx=newx[, attr(object[[i]], "preSelectVars"), drop=FALSE])
        cat("Ensemble", round(i/ensembleSize, 4)*100, "%", "\n")
    }
    else{
      predMat[, i] <- predict(object=object[[i]], newx=newx)
      cat("Ensemble", round(i/ensembleSize, 4)*100, "%", "\n")
    }
  }
  return(predMat)
}

####################################
# Ensemble predictions saved on disk

predict.KDSNensembleDisk <- function (object, newx, ...) {
  ensembleSize <- length(object)
  predMat <- matrix(NA, nrow=nrow(newx), ncol=ensembleSize)
  ensembleModel <- NULL
  for(i in 1:ensembleSize) {
    
    # Load ensemble model
    load(object[i])
    
    varPreSelectInput <- any(names(attributes(ensembleModel))=="preSelectVars")
    if(varPreSelectInput) {
      predMat[, i] <- predict(object=ensembleModel, newx=newx) 
#                              newx=newx[, attr(ensembleModel, "preSelectVars"), drop=FALSE])
      cat("Ensemble", round(i/ensembleSize, 4)*100, "%", "\n")
    }
    else{
      predMat[, i] <- predict(object=ensembleModel, newx=newx)
      cat("Ensemble", round(i/ensembleSize, 4)*100, "%", "\n")
    }
    
    # Remove model from workspace
    rm(ensembleModel)
    
  }
  return(predMat)
}
