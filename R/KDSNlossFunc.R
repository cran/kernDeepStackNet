##################################################
# cross validation loss

devStandard <- function (preds, ytest, RMSE=TRUE) {
  if(RMSE){
    return(sqrt(sum( (ytest - preds)^2 ) / length(ytest) ))
  }
  else{
    return(sum((ytest-preds)^2))
  }
}

#####################
# Main loss function
#####################

lossCvKDSN <- function (parOpt, y, X, levels, cvIndex, seedW=NULL, lossFunc=devStandard,
                        varSelect=rep(FALSE, levels), varRanking=NULL, 
                        alpha=rep(0, levels), dropHidden=rep(FALSE, levels),
                        seedDrop=NULL, standX=TRUE, standY=FALSE,
                        namesParOpt=rep(c("dim", "sigma", "lambdaRel"), levels)) {
  # Checks
  lenParOpt <- length(parOpt)
  stopifnot((lenParOpt %% levels) ==0)
  stopifnot(length(seedW)==levels)
  stopifnot(length(dropHidden)==levels)
  stopifnot(length(alpha)==levels)
  stopifnot(length(varSelect)==levels)

#   # Extract parameters
#   varExFreq <- parOpt [seq(1, lenParOpt, 5)]
#   Dim <- round (parOpt [seq(2, lenParOpt, 5)])
#   sigma <- parOpt [seq(3, lenParOpt, 5)]
#   dropHiddenProb <- parOpt [seq(4, lenParOpt, 5)]
#   lambdaRel <- parOpt [seq(5, lenParOpt, 5)]
  
  # Extract parameters
  varExFreq <- parOpt[grep("select", namesParOpt)]
  Dim <- round(parOpt[grep("dim", namesParOpt)])
  sigma <- parOpt[grep("sigma", namesParOpt)]
  dropHiddenProb <- parOpt[grep("dropProb", namesParOpt)]
  lambdaRel <- parOpt[grep("lambdaRel", namesParOpt)]
  # Enlarge with standard parameters to use correct values in fitting
  if(sum(varSelect)!=levels) {
    tempVarExFreq <- vector("numeric", levels)
    tempVarExFreq[varSelect] <- varExFreq
    tempVarExFreq[!varSelect] <- NA
    varExFreq <- tempVarExFreq
  }
  if(sum(dropHidden)!=levels) {
    tempDropHiddenProb <- vector("numeric", levels)
    tempDropHiddenProb[dropHidden] <- dropHiddenProb
    tempDropHiddenProb[!dropHidden] <- NA
    dropHiddenProb <- tempDropHiddenProb
  }
  
  # Compute deviance on test data
  cvSamples <- length(cvIndex)
  vecLoss <- vector("numeric", cvSamples)
  for(j in 1:cvSamples) {
    
    # Fit model
    fitKDSN <- fitKDSN (y=y [cvIndex [[j]] ], 
                        X=X [cvIndex [[j]], , drop=FALSE], 
                        levels=levels, Dim=Dim, sigma=sigma, lambdaRel=lambdaRel, 
                        alpha=alpha, info=FALSE, seedW=seedW, standX=standX, standY=standY,
                        varSelect=varSelect, varRanking=varRanking, varExFreq=varExFreq,
                        dropHidden=dropHidden, dropHiddenProb=dropHiddenProb, 
                        seedDrop=seedDrop)
    
    # Compute predicted values
    predVal <- predict(fitKDSN, newx=X [-cvIndex [[j]], , drop=FALSE])

    # Compute loss
    vecLoss [j] <- lossFunc(preds=predVal, ytest=y [-cvIndex [[j]] ])
  }
  
  # Calculate average loss and output
  Loss <- mean(vecLoss)
  attributes(Loss) <- list(Fit=fitKDSN)
  return(Loss)
}

# Help function
# Computes the predictive log probability for model based optimization
predLogProb <- function(predMean, predSigma, y, X) {
  sum(-log(predSigma)/2 - 
        (y-predMean)^2 / (2 * predSigma) - 
        log(2*pi) / 2)
}

#####################################
# Loss function for shared parameters
#####################################

lossSharedCvKDSN <- function (parOpt, y, X, cvIndex, seedW=NULL, lossFunc=devStandard,
                        varSelectShared=FALSE, varRanking=NULL, 
                        alphaShared=0, dropHiddenShared=FALSE,
                        seedDrop=NULL, standX=TRUE, standY=FALSE,
                        namesParOpt=c("levels", "dim", "sigma", "lambdaRel")) {

  # Extract tuning parameters
  levels <- round(parOpt[grep("levels", namesParOpt)])
  if(varSelectShared) {
    varExFreq <- rep(parOpt[grep("select", namesParOpt)], levels)
  }
  else{
    varExFreq <- rep(NA, levels)
  }
  Dim <- rep(round(parOpt[grep("dim", namesParOpt)]), levels)
  sigma <- rep(parOpt[grep("sigma", namesParOpt)], levels)
  if(dropHiddenShared) {
    dropHiddenProb <- rep(parOpt[grep("dropProb", namesParOpt)], levels)
  }
  else{
    dropHiddenProb <- rep(NA, levels)
  }
  lambdaRel <- rep(parOpt[grep("lambdaRel", namesParOpt)], levels)

  # Fixed parameters
  alpha <- rep(alphaShared, levels)
  dropHidden <- rep(dropHiddenShared, levels)
  varSelect <- rep(varSelectShared, levels)
  
  # Compute deviance on test data
  cvSamples <- length(cvIndex)
  vecLoss <- vector("numeric", cvSamples)
  for(j in 1:cvSamples) {
    
    # Fit model
    fitKDSN <- fitKDSN (y=y [cvIndex [[j]] ], 
                        X=X [cvIndex [[j]], , drop=FALSE], 
                        levels=levels, Dim=Dim, sigma=sigma, lambdaRel=lambdaRel, 
                        alpha=alpha, info=FALSE, seedW=seedW, standX=standX, standY=standY,
                        varSelect=varSelect, varRanking=varRanking, 
                        varExFreq=rep(varExFreq, levels),
                        dropHidden=dropHidden, 
                        dropHiddenProb=dropHiddenProb, 
                        seedDrop=seedDrop)
    
    # Compute predicted values
    predVal <- predict(fitKDSN, newx=X [-cvIndex [[j]], , drop=FALSE])
    
    # Compute loss
    vecLoss [j] <- lossFunc(preds=predVal, ytest=y [-cvIndex [[j]] ])
  }
  
  # Calculate average loss and output
  Loss <- mean(vecLoss)
  attributes(Loss) <- list(Fit=fitKDSN)
  return(Loss)
}

# Evaluate loss function on external test data
lossSharedTestKDSN <- function (parOpt, y, X, yTest, Xtest, 
                                seedW=NULL, lossFunc=devStandard,
                              varSelectShared=FALSE, varRanking=NULL, 
                              alphaShared=0, dropHiddenShared=FALSE,
                              seedDrop=NULL, standX=TRUE, standY=FALSE,
                              namesParOpt=c("levels", "dim", "sigma", "lambdaRel")) {
  
  # Extract tuning parameters
  levels <- round(parOpt[grep("levels", namesParOpt)])
  if(varSelectShared) {
    varExFreq <- rep(parOpt[grep("select", namesParOpt)], levels)
  }
  else{
    varExFreq <- rep(NA, levels)
  }
  Dim <- rep(round(parOpt[grep("dim", namesParOpt)]), levels)
  sigma <- rep(parOpt[grep("sigma", namesParOpt)], levels)
  if(dropHiddenShared) {
    dropHiddenProb <- rep(parOpt[grep("dropProb", namesParOpt)], levels)
  }
  else{
    dropHiddenProb <- rep(NA, levels)
  }
  lambdaRel <- rep(parOpt[grep("lambdaRel", namesParOpt)], levels)
  
  # Fixed parameters
  alpha <- rep(alphaShared, levels)
  dropHidden <- rep(dropHiddenShared, levels)
  varSelect <- rep(varSelectShared, levels)
  
  # Compute deviance on test data

  # Fit model
  fitKDSN <- fitKDSN (y=y, X=X, levels=levels, Dim=Dim, sigma=sigma, 
                      lambdaRel=lambdaRel, alpha=alpha, info=FALSE, 
                      seedW=seedW, standX=standX, standY=standY,
                      varSelect=varSelect, varRanking=varRanking, 
                      varExFreq=rep(varExFreq, levels),
                      dropHidden=dropHidden, 
                      dropHiddenProb=dropHiddenProb, 
                      seedDrop=seedDrop)
  
  # Compute predicted values
  predVal <- predict(fitKDSN, newx=Xtest[, , drop=FALSE])
  
  # Compute loss loss and output
  Loss <- lossFunc(preds=predVal, ytest=yTest)
  attributes(Loss) <- list(Fit=fitKDSN)
  return(Loss)
}

####################################
# Loss with GCV approximation
# Reimplemented from old source code
####################################

lossApprox <- function (parOpt, y, X, levels, seedW=NULL,
                        varSelect=rep(FALSE, levels), varRanking=NULL, 
                        alpha=rep(0, levels), dropHidden=rep(FALSE, levels),
                        seedDrop=NULL, standX=TRUE, standY=FALSE,
                        namesParOpt=rep(c("dim", "sigma", "lambdaRel"), levels),
                        gammaPar=1) {
  # Checks
  lenParOpt <- length(parOpt)
  stopifnot((lenParOpt %% levels) ==0)
  stopifnot(length(seedW)==levels)
  stopifnot(length(dropHidden)==levels)
  stopifnot(length(alpha)==levels)
  stopifnot(length(varSelect)==levels)

  # Extract parameters
  varExFreq <- parOpt[grep("select", namesParOpt)]
  Dim <- round(parOpt[grep("dim", namesParOpt)])
  sigma <- parOpt[grep("sigma", namesParOpt)]
  dropHiddenProb <- parOpt[grep("dropProb", namesParOpt)]
  lambdaRel <- parOpt[grep("lambdaRel", namesParOpt)]
  # Enlarge with standard parameters to use correct values in fitting
  if(sum(varSelect)!=levels) {
    tempVarExFreq <- vector("numeric", levels)
    tempVarExFreq[varSelect] <- varExFreq
    tempVarExFreq[!varSelect] <- NA
    varExFreq <- tempVarExFreq
  }
  if(sum(dropHidden)!=levels) {
    tempDropHiddenProb <- vector("numeric", levels)
    tempDropHiddenProb[dropHidden] <- dropHiddenProb
    tempDropHiddenProb[!dropHidden] <- NA
    dropHiddenProb <- tempDropHiddenProb
  }
  
  # Compute deviance on test data
  # Fit model
  KDSNfit <- fitKDSN (y=y, X=X, levels=levels, Dim=Dim, sigma=sigma, lambdaRel=lambdaRel, 
                      alpha=alpha, info=FALSE, seedW=seedW, standX=standX, standY=standY,
                      varSelect=varSelect, varRanking=varRanking, varExFreq=varExFreq,
                      dropHidden=dropHidden, dropHiddenProb=dropHiddenProb, 
                      seedDrop=seedDrop)
    
  # Compute predicted values
  fittedVal <- predict(KDSNfit, newx=X)

  # Calculate weight matrix of observations
  varFunc <- varMu (mu=fittedVal)
  linkDeriv <- gDerivMu (mu=fittedVal)
  Wdiag <- calcWdiag (varMu=varFunc, gDerivMu=linkDeriv)
  
  # Extract deviance of last level
  fitDev <- devStandard(preds=fittedVal, ytest=c(y))
  
  # Calculate influence (hat) matrix
  trAmat <- calcTrAFast (X=t(KDSNfit$Output$rftX [[KDSNfit$Input$levels]] $Z), 
                         w=Wdiag, lambda=KDSNfit$Output$lambdaActual[levels])
  
  # Compute loss
  Loss <- lossGCV (n=nrow(X), Dev=fitDev, trA=trAmat, gammaPar=gammaPar)
  
  # Numerical check of loss
  CheckCondition <- is.infinite(Loss) | 
    is.null(Loss) | 
    is.nan(Loss) |
    is.na(Loss)
  if(any(CheckCondition)) {
    Loss[CheckCondition] <- max(Loss, na.rm=TRUE)*2
  }
  
  # Output
  attributes(Loss) <- list(Fit=KDSNfit)
  return(Loss)
}

# GCV loss
lossGCV <- function (n, Dev, trA, gammaPar=1) {
  stopifnot((n - gammaPar * trA) >= 0)
  n * Dev / (n - gammaPar * trA)^2
}

# Calculate W matrix in P-IRLS
calcWdiag <- function (varMu, gDerivMu) {
  c((varMu * gDerivMu^2) ^ (-1/2))
}

# Variance function of glm
varMu <- function (mu) {
  return(rep(1, length(mu)))
}

# Derivative of link function of glm 
gDerivMu <- function (mu) {
  return(rep(1, length(mu)))
}

calcTrA <- function (W, X, lambda=0) {
  
  # Matrix calculation
  B_mat <- W %*% X
  
  # QR decomposition
  R_mat <- qr.R(qr(B_mat)) # p x p
  E_mat <- diag(rep(sqrt(lambda), dim(X) [2])) # p x p
  R_E <- rbind(R_mat, E_mat)
  
  # Singular value decomposition
  U_mat_1 <- svd(R_E)$u [dim(R_mat) [1], ]
  
  # Compute trace of hat matrix A
  trA <- sum(diag(U_mat_1 %*% t(U_mat_1)))
  return(trA)
}

# Fast trace computation
calcTrAFast <- function(X, w, lambda=0) {
  Xtilde <- X * w
  XtildeT <- crossprodRcpp(Xtilde)[[1]]
  eigenVal <- getEigenValuesRcpp(XtildeT)
  sum(eigenVal/(eigenVal + lambda))
}
