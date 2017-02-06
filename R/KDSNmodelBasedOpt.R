####################################################
# Model based optimization based on Kriging with EGO

mbo1d <- function (model, fun, nsteps, lower, upper, parinit, isoInput, maxRunsMult=1, 
                   repMult=1, tol_input=.Machine$double.eps^0.25, addInfo=TRUE, 
                   # nCores=1, 
                   envir=parent.frame(), EIopt="1Dmulti", GenSAmaxCall=100, timeAlloc="constant",
                   EItype="EI") {

  # Check to increase bounds
  checkFuncBody <- as.character(functionBody(fun))[1]
  # if( (checkFuncBody=="lossSharedCvKDSN" | 
  #     checkFuncBody=="lossSharedTestKDSN") & (lower[1]==upper[1]) ) {
  #   upper[1] <- upper[1] + 0.1
  # }
  for(i in 1:nsteps) {
    
    # Time allocation
    # old code
    # function (x) -EImod(x=x, model=model)
    if(timeAlloc=="constant") {
      newNoiseVar <- model@covariance@nugget / (nsteps-i+1)
    }
    if(timeAlloc=="zero") {
      newNoiseVar <- 0
    }
    
    if(EIopt=="1Dmulti") {
      
      if( (checkFuncBody=="lossSharedCvKDSN" | 
           checkFuncBody=="lossSharedTestKDSN") & (lower[1]==1 & upper[1]==1) ) {
        # Expected quantile improvement criterion
        if(EItype=="EQI") {
          opt_res <- optimize1dMulti (f_input=function (x) -EQI(x=x, model=model, 
                                                                new.noise.var = newNoiseVar, beta=0.9),
                                      lower=lower[-1], upper=upper[-1],
                                      maxRuns=length(parinit[-1])*maxRunsMult, 
                                      repetitions=length(parinit[-1])*repMult, 
                                      tol_input=tol_input, x_0=parinit[-1], addInfo=FALSE, 
                                      nCores=1, envir=parent.frame(), directUse=FALSE, OptTypePar="mbo1d")
        }
        # Expected improvement criterion
        if(EItype=="EI") {
          opt_res <- optimize1dMulti (f_input=function (x) -EI(x=x, model=model),
                                      lower=lower[-1], upper=upper[-1],
                                      maxRuns=length(parinit[-1])*maxRunsMult, 
                                      repetitions=length(parinit[-1])*repMult, 
                                      tol_input=tol_input, x_0=parinit[-1], addInfo=FALSE, 
                                      nCores=1, envir=parent.frame(), directUse=FALSE, OptTypePar="mbo1d")
        }
      
        # Construct new design
        new_func_val <- c(model@y, fun(c(1, opt_res$minimum)))
        new_design <- cbind(1, rbind(model@X, opt_res$minimum))
      }
      else{
        # Optimize improvement criterion
        if(EItype=="EQI") {
          opt_res <- optimize1dMulti (f_input=function (x) -EQI(x=x, model=model, 
                                                                new.noise.var = newNoiseVar, beta=0.9),
                                      lower=lower, upper=upper,
                                      maxRuns=length(parinit)*maxRunsMult, repetitions=length(parinit)*repMult, 
                                      tol_input=tol_input, x_0=parinit, addInfo=FALSE, 
                                      nCores=1, envir=parent.frame(), directUse=FALSE, OptTypePar="mbo1d")
        }
        if(EItype=="EI") {
          opt_res <- optimize1dMulti (f_input=function (x) -EI(x=x, model=model),
                                      lower=lower, upper=upper,
                                      maxRuns=length(parinit)*maxRunsMult, repetitions=length(parinit)*repMult, 
                                      tol_input=tol_input, x_0=parinit, addInfo=FALSE, 
                                      nCores=1, envir=parent.frame(), directUse=FALSE, OptTypePar="mbo1d")
        }
      
        # Construct new design
        new_func_val <- c(model@y, fun(opt_res$minimum))
        new_design <- rbind(model@X, opt_res$minimum)
      }
    }
    
    if(EIopt=="GenSA") {
      if( (checkFuncBody=="lossSharedCvKDSN" | 
           checkFuncBody=="lossSharedTestKDSN") & (lower[1]==1 & upper[1]==1) ) {
        # Optimize expected improvement criterion
        if(EItype=="EQI") { 
          opt_res <- GenSA(par=parinit[-1], fn=function (x) -EQI(x=x, model=model, 
                                                                 new.noise.var = newNoiseVar, beta=0.9), 
                                             lower=lower[-1], 
                                             upper=upper[-1], 
                                             control=list(max.call=GenSAmaxCall*length(parinit[-1])))
          }
        
        if(EItype=="EI") {
          opt_res <- GenSA(par=parinit[-1], fn=function (x) -EI(x=x, model=model), 
                           lower=lower[-1], 
                           upper=upper[-1], 
                           control=list(max.call=GenSAmaxCall*length(parinit[-1])))
        }

        # Construct new design
        new_func_val <- c(model@y, fun(c(1, opt_res$par)))
        new_design <- cbind(1, rbind(model@X, opt_res$par))
      }
      else{
        # Optimize expected improvement criterion
        if(EItype=="EQI") {
          opt_res <- GenSA(par=parinit, fn=function (x) -EQI(x=x, model=model, 
                                                             new.noise.var = newNoiseVar, beta=0.9), 
                           lower=lower, 
                           upper=upper, 
                           control=list(max.call=GenSAmaxCall*length(parinit)))
        }
        
        if(EItype=="EI") {
          opt_res <- GenSA(par=parinit, fn=function (x) -EI(x=x, model=model), 
                           lower=lower, 
                           upper=upper, 
                           control=list(max.call=GenSAmaxCall*length(parinit)))
        }

        # Construct new design
        new_func_val <- c(model@y, fun(opt_res$par))
        new_design <- rbind(model@X, opt_res$par)
      }
    }

    # Reevaluate Kriging model
    parinit <- new_design [which.min(new_func_val), ]
    # Correction for setting maxLevel=1
    if((checkFuncBody=="lossSharedCvKDSN" | 
        checkFuncBody=="lossSharedTestKDSN") & 
       var(new_design[, 1])==0) {
      model <- km(design=new_design[, -1], response=new_func_val, control=list(trace=FALSE), nugget.estim=TRUE, iso=isoInput)
    }
    else{
      model <- km(design=new_design, response=new_func_val, control=list(trace=FALSE), nugget.estim=TRUE, iso=isoInput)
    }
    if(addInfo) {cat("N_step", i, "of", nsteps, "\n")}
  }

  # Return best solution
  Index <- which.min(new_func_val)
  result <- new_design [Index, ]
  attributes(result) <- list(funcVal=new_func_val [Index], MBOprogress=cbind(Loss=new_func_val, new_design))
  return(result)
}

# Guidelines:
# Use default matern 5/2 correlation structure to avoid numerical problems during invertation of matrices
# Use number of steps and starting points proportional to dimension D, e. g. nsteps= startPoints = D * 10

# Input Arguments
# loss_func: Loss function with only a multivariate vector as input and a scalar output value!
# n_steps: Number of initial points and number of steps of EGO
# lower_bounds: Vector of lower bounds of variables
# upper_bounds: Vector of upper bounds of variables
# x_start: Starting value for optimization
# Output list with names
# par: Best parameters
# value: Function value at the best parameters

mboAll <- function (loss_func, n_steps, initDesign, lower_bounds, upper_bounds, x_start, isoInput=FALSE, addInfo=TRUE,
                    maxRunsMult=1, repMult=1, tol_input=.Machine$double.eps^0.25, envir=parent.frame(),
                    # nCores=1, 
                    metaModelSelect=TRUE, EIopt="1Dmulti", GenSAmaxCall=100, timeAlloc="constant",
                    EItype="EI") {
#  if(nCores<1){stop("Please specify an integer number greater or equal to one as the number of threads!")}
  # Generate design and transform 
  dimVec <- length(x_start)
  LHS_design <- maximinLHS(n=initDesign, k=dimVec)
  LHS_design <- sapply(1:dimVec, function (x) (upper_bounds[x] - lower_bounds[x]) * LHS_design [, x] + lower_bounds[x])
  if(addInfo) {cat("LHS_design", "\n")}
  
  # Evaluate loss function on design
  func_start <- vector("numeric", initDesign)
  localEnvir <- environment()
#  if(nCores==1) {
    for(i in 1:initDesign) {
      func_start[i] <- loss_func(LHS_design[i, ])
      if(addInfo) {cat("Iter", i, "of", initDesign, "\n")}
    }
#  }
  # else {
  #   cl <- makeCluster(nCores)
  #   clusterExport(cl = cl, varlist=c("loss_func", "LHS_design", "initDesign"), envir = localEnvir)
  #   # Necessary to access variables in function
  #   tempVec <- x_start
  #   tempVec[seq(1, length(tempVec), 5)] <- 1
  #   loss_func(tempVec)[1]
  #   func_start <- parSapplyLB(cl=cl, X=1:initDesign, 
  #                             FUN=function(x) loss_func(LHS_design[x, ])[1])
  #   stopCluster(cl=cl)
  # }
  if(addInfo) {cat("LHS_design func eval", "\n")}
  
  checkFuncBody <- as.character(functionBody(loss_func))[1]
  if(metaModelSelect) {
  # Metamodel validation and Kriging estimation
  # Trend: constant
  # covtype: "gauss", "matern5_2", "matern3_2", "exp" or "powexp"
  covtypeInput <- c("gauss", "matern5_2", "matern3_2", "exp", "powexp")
  perfMeasure <- vector("numeric", length=5)
  km_select <- vector("list", 5)
  #  if(nCores==1) {
    for(j in 1:5) {
      # Correction for setting maxLevel=1
      if((checkFuncBody=="lossSharedCvKDSN" | 
          checkFuncBody=="lossSharedTestKDSN") & 
         var(LHS_design[, 1])==0) {
        km_select[[j]] <- km(covtype=covtypeInput[j], design=LHS_design[, -1], response=func_start, control=list(trace=FALSE), nugget.estim=TRUE, iso=isoInput)
      }
      else{
        km_select[[j]] <- km(covtype=covtypeInput[j], design=LHS_design, response=func_start, control=list(trace=FALSE), nugget.estim=TRUE, iso=isoInput)
      }
      LOOmuSigma <- leaveOneOut.km(model=km_select[[j]], type="UK", trend.reestim=TRUE)
      perfMeasure [j] <- predLogProb (predMean=LOOmuSigma$mean, predSigma=LOOmuSigma$sd^2, y=func_start, X=LHS_design)
    }
    km_base <- km_select [[which.max(perfMeasure)]]
#  }
  # else {
  #   tempFunction <- function(j) {
  #     km_select1 <- km(covtype=covtypeInput[j], design=LHS_design, response=func_start, control=list(trace=FALSE), nugget.estim=TRUE, iso=isoInput)
  #     LOOmuSigma <- leaveOneOut.km(model=km_select1, type="UK", trend.reestim=TRUE)
  #     perfMeasure1 <- predLogProb (predMean=LOOmuSigma$mean, predSigma=LOOmuSigma$sd^2, y=func_start, X=LHS_design)
  #     return(list(Perf=perfMeasure1, Model=km_select1))
  #   }
  #   cl <- makeCluster(ifelse(nCores<=5, nCores, 5))
  #   clusterExport(cl = cl, varlist=c("tempFunction", "func_start", "isoInput", "covtypeInput"), envir = localEnvir)
  #   tempResults <- parLapply(cl=cl, X=1:5, fun=tempFunction)
  #   stopCluster(cl=cl)
  #   maxInd <- which.max(sapply(1:5, function(x) tempResults[[x]]$Perf))
  #   km_base <- tempResults [[maxInd]]$Model
  # }
    if(addInfo) {cat("Kriging model selection", "\n")}
  }
  else{
    # Correction for setting maxLevel=1
    if((checkFuncBody=="lossSharedCvKDSN" | 
        checkFuncBody=="lossSharedTestKDSN") & 
       var(LHS_design[, 1])==0) {
      km_base <- km(design=LHS_design[, -1], response=func_start, control=list(trace=FALSE), nugget.estim=TRUE, iso=isoInput)
    }
    else{
      km_base <- km(design=LHS_design, response=func_start, control=list(trace=FALSE), nugget.estim=TRUE, iso=isoInput)
    }
  }

  # Optimize
  ego_result <- mbo1d(model=km_base, fun=loss_func, nsteps=n_steps, lower=lower_bounds, upper=upper_bounds, parinit=x_start, 
                      isoInput=isoInput, addInfo=addInfo,
                      maxRunsMult=maxRunsMult, repMult=repMult, tol_input=tol_input, 
                      EIopt=EIopt, GenSAmaxCall=GenSAmaxCall, timeAlloc=timeAlloc,
                      EItype=EItype)
  if(addInfo){cat("MBO done", "\n")}
  # Output
  return(list(par=c(ego_result), value=attr(ego_result, "funcVal"), 
              MBOprogress=attr(ego_result, "MBOprogress")))
}

###############################
# Main tuning function of KDSN
###############################

tuneMboLevelCvKDSN <- function (y, X, levels=1, alpha=rep(0, levels), fineTuneIt=100, 
                                nStepMult=20, designMult=10, 
                                dimMax=round(sqrt(dim(X)[1])/2), addInfo=TRUE, # nCores=1,
                                maxRunsMult=1, repMult=1, tol_input=.Machine$double.eps^0.25, 
                                cvIndex, lossFunc=devStandard, EIopt="1Dmulti", GenSAmaxCall=100,
                                varSelect=rep(FALSE, levels), rdcRep=1, dropHidden=rep(FALSE, levels),
                                standX=TRUE, standY=FALSE, timeAlloc="constant", varPreSelect=FALSE,
                                varPreSelpopSize=100, varPreSelMaxiter=100, EItype="EQI") {

  # if(nCores<1){stop("Please specify an integer number greater or equal to one as the number of threads!")}
  #   if(nCores>1) {
  #     # Initialize cluster
  #     cl <- makeCluster(nCores)
  #   }
  
  preSelectIndices <- NULL
  if(varPreSelect) {
    # Exclude columns prior variable selection with zero variance
    zeroSdVars <- which(colSds(X)==0)
    convertIndices <- setdiff(1:ncol(X), zeroSdVars)
    if(length(zeroSdVars)!=0) {
      X <- X[, -zeroSdVars, drop=FALSE]
    }
    
    # RDC variable selection based on all available covariates
    preSelectIndices <- rdcVarSelSubset(x=X, y=y, k=20, 
                                        s=1/6, f=sin, seedX=1:rdcRep, seedY=-c(1:rdcRep), 
                                        rdcRep=rdcRep, popSize=varPreSelpopSize, 
                                        maxiter=varPreSelMaxiter)

    # Select important variables
    X <- X[, preSelectIndices, drop=FALSE]
    
    # Convert indices to original scale for future predictions
    preSelectIndices <- convertIndices[preSelectIndices]
    
    if(addInfo) {cat("RDC variable pre-selection", "\n")}
  }
  
  rdcIndizes <- NULL
  if(any(varSelect)) {
      rdcIndizes <- rdcVarOrder(x=X, y=y, cutoff=0, # nCores=nCores, 
                              seedX=1:rdcRep, seedY=-(1:rdcRep), rdcRep=rdcRep)
      if(addInfo) {cat("RDC marginal variable ordering", "\n")}
  }
  
  # Initialize parameters
  n <- dim(X)[1]
  
  # Initialize starting vector of hyperparameters
  varExFreqStart <- 0.5
  MaxMatEntries <- .Machine$integer.max
  dimStart <- round ((dimMax+1) / 2)
  if((dimMax*2*n) > MaxMatEntries) {
    dimMax <- floor(MaxMatEntries/n/2)
    dimStart <- round (sqrt(dimMax*2)/2)
  }
  quantEuklid <- quantile(c(dist(robustStandard(X))^2), probs=c(0, 0.5, 1))
  # Correct for possible zero as lower bound of sigma
  if(any(quantEuklid==0)) {
    quantEuklid[quantEuklid==0] <- .Machine$double.eps^0.5
  }
  sigmaStart <- quantEuklid [2]
  lambdaRelStart <- 0.5
  dropHiddenProbStart <- 0.5

# old code
#   x_start <- c(varExFreqStart, dimStart, sigmaStart, dropHiddenProbStart, lambdaRelStart)
#   x_new <- rep(x_start, levels)
#
#   # Initialize bounds of hyperparameters
#   interval_matrix_start <- matrix(c(0, 1, 1, dimMax, quantEuklid [1], quantEuklid [3], 0.025, 0.975, 0, 0.99), nrow=2, ncol=5)
#   interval_matrix <- interval_matrix_start [, rep(1:dim(interval_matrix_start)[2], levels)]

  # Consider four different cases:
  # 1: Variable selection and dropout
  # 2: Dropout
  # 3: Variable selection
  # 4: Neither variable selection nor dropout
  
  tempInd <- as.numeric(varSelect) + as.numeric(dropHidden) + 3
  noPar <- sum(tempInd)
  cumIndex <- c(0, cumsum(tempInd))
  Indices <- lapply(2:length(cumIndex), function(x) (1+cumIndex[x-1]):(cumIndex[x]))
  x_new <- vector("numeric", noPar)
  lower_bounds <- vector("numeric", noPar)
  upper_bounds <- vector("numeric", noPar)
  nameVec <- vector("character", noPar)
  for(l in 1:levels) {
    if(varSelect[l] & dropHidden[l]){
      x_new[Indices[[l]] ] <- c(varExFreqStart, dimStart, sigmaStart, dropHiddenProbStart, lambdaRelStart)
      lower_bounds[Indices[[l]] ] <- c(0, 1, quantEuklid [1], 0.025, 0)
      upper_bounds[Indices[[l]] ] <- c(1, dimMax, quantEuklid [3], 0.975, 0.99)
      nameVec[Indices[[l]] ] <- c("select", "dim", "sigma", "dropProb", "lambdaRel")
    }
    if(!varSelect[l] & dropHidden[l]){
      x_new[Indices[[l]] ] <- c(dimStart, sigmaStart, dropHiddenProbStart, lambdaRelStart)
      lower_bounds[Indices[[l]] ] <- c(1, quantEuklid [1], 0.025, 0)
      upper_bounds[Indices[[l]] ] <- c(dimMax, quantEuklid [3], 0.975, 0.99)
      nameVec[Indices[[l]] ] <- c("dim", "sigma", "dropProb", "lambdaRel")
    }
    if(varSelect[l] & !dropHidden[l]){
      x_new[Indices[[l]] ] <- c(varExFreqStart, dimStart, sigmaStart, lambdaRelStart)
      lower_bounds[Indices[[l]] ] <- c(0, 1, quantEuklid [1], 0)
      upper_bounds[Indices[[l]] ] <- c(1, dimMax, quantEuklid [3], 0.99)
      nameVec[Indices[[l]] ] <- c("select", "dim", "sigma", "lambdaRel")
    }
    if(!varSelect[l] & !dropHidden[l]){
      x_new[Indices[[l]] ] <- c(dimStart, sigmaStart, lambdaRelStart)
      lower_bounds[Indices[[l]] ] <- c(1, quantEuklid [1], 0)
      upper_bounds[Indices[[l]] ] <- c(dimMax, quantEuklid [3], 0.99)
      nameVec[Indices[[l]] ] <- c("dim", "sigma", "lambdaRel")
    }
  }

  # Tune
  optVal <- mboAll (loss_func=function (x) lossCvKDSN (parOpt=x, y=y, X=X, cvIndex=cvIndex, levels=levels,  
                                                     seedW=seq(0, (levels-1), 1), lossFunc=lossFunc,
                                                     varSelect=varSelect, varRanking=rdcIndizes, alpha=alpha,
                                                     dropHidden=dropHidden, seedDrop=seq(0, -(levels-1), -1), 
                                                     standX=standX, standY=standY, namesParOpt=nameVec), 
                    n_steps=(length(x_new)+3)*nStepMult, initDesign=(length(x_new)+3)*designMult, 
                    lower_bounds=lower_bounds, upper_bounds=upper_bounds, x_start=x_new,
                    # nCores=nCores, 
                    addInfo=addInfo, 
                    maxRunsMult=maxRunsMult, repMult=repMult, tol_input=tol_input, 
                    EIopt=EIopt, GenSAmaxCall=GenSAmaxCall, timeAlloc=timeAlloc,
                    EItype=EItype)
  x_new <- optVal$par
  
  # Fine tune random fourier transformation weights
  # Reproduceability is ensured with seed generation
  if(fineTuneIt > 0) {
    fineTune <- vector("numeric", fineTuneIt)
    seedGuess <- matrix(sample.int(.Machine$integer.max, size=fineTuneIt * levels) * sample(c(-1, 1), size=fineTuneIt * levels, replace=TRUE), 
                        nrow=fineTuneIt, ncol=levels)
    seedDropGuess <- matrix(sample.int(.Machine$integer.max, size=fineTuneIt * levels) * sample(c(-1, 1), size=fineTuneIt * levels, replace=TRUE), 
                        nrow=fineTuneIt, ncol=levels)
    localEnvir <- environment()
#    if(nCores==1) {
      for(i in 1:fineTuneIt) {
        fineTune[i] <- lossCvKDSN(parOpt=x_new, y=y, X=X, cvIndex=cvIndex, 
                                  seedW=seedGuess [i, ], lossFunc=lossFunc, 
                                  varSelect=varSelect, varRanking=rdcIndizes, 
                                  alpha=alpha, levels=levels,
                                  dropHidden=dropHidden,
                                  seedDrop=seedDropGuess[i, ], 
                                  standX=standX, standY=standY, namesParOpt=nameVec)[1]
        if(addInfo) {cat("tuneKDSN", "FineTune =", i, "\n")}
      }
#    }
    # else {
    #   cl <- makeCluster(nCores)
    #   clusterExport(cl = cl, varlist=c("lossCvKDSN", "x_new", "y", "X", "seedGuess", "fineTuneIt"), 
    #                 envir = localEnvir)
    #   fineTune <- parSapply(cl=cl, X=1:fineTuneIt, 
    #                         FUN=function(i) lossCvKDSN(parOpt=x_new, y=y, X=X, cvIndex=cvIndex,  
    #                                                    seedW=seedGuess [i, ], lossFunc=lossFunc,
    #                                                    varSelect=varSelect, varRanking=rdcIndizes,
    #                                                    alpha=alpha, levels=levels,
    #                                                    dropHidden=dropHidden,
    #                                                    seedDrop=seedDropGuess[i, ],
    #                                                    standX=standX, standY=standY, namesParOpt=nameVec)[1])
    #   stopCluster(cl=cl)
    #   if(addInfo) {cat("FineTuning", "\n")}
    # }
    minIndex <- which.min(fineTune)
  }
  else {
    fineTune <- Inf
    minIndex <- 1
  }
  
  # Output
  # Refit best model
  varExFreq <- x_new[grep("select", nameVec)]
  Dim <- round(x_new[grep("dim", nameVec)])
  sigma <- x_new[grep("sigma", nameVec)]
  dropHiddenProb <- x_new[grep("dropProb", nameVec)]
  lambdaRel <- x_new[grep("lambdaRel", nameVec)]
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

  if(fineTune[minIndex] < optVal$value){
    finalModel <- fitKDSN(y = y, X = X, levels = levels, Dim = Dim, 
                          sigma = sigma, lambdaRel = lambdaRel, alpha = alpha, 
                          info = FALSE, seedW = seedGuess [minIndex, ], 
                          standX = standX, standY=standY,
                          varSelect=varSelect, varRanking=rdcIndizes,
                          varExFreq=varExFreq, 
                          dropHidden=dropHidden, 
                          dropHiddenProb=dropHiddenProb,
                          seedDrop=seedDropGuess[minIndex, ])
    attr(finalModel, which="Loss") <- fineTune[minIndex]
  }
  else{
    finalModel <- fitKDSN(y = y, X = X, levels = levels, Dim = Dim, 
                          sigma = sigma, lambdaRel = lambdaRel, alpha = alpha, 
                          info = FALSE, seedW = seq(0, (levels-1), 1), 
                          standX = standX, standY=standY,
                          varSelect=varSelect, varRanking=rdcIndizes,
                          varExFreq=varExFreq,
                          dropHidden=dropHidden, 
                          dropHiddenProb=dropHiddenProb,
                          seedDrop=seq(0, -(levels-1), -1))
    attr(finalModel, which="Loss") <- optVal$value
  }
  # Include Loss score, loss function, cvIndex and preselected variables as attributes
  attr(finalModel, which="LossFunc") <- lossFunc
  attr(finalModel, which="cvIndex") <- cvIndex
  if(varPreSelect) {
    finalModel$Output$varPreSelect <- TRUE
    attr(finalModel, which="preSelectVars") <- preSelectIndices
  }
  # Include progress of MBO optimization
  dimnames(optVal$MBOprogress)[[2]] <- c("Loss", nameVec)
  attr(finalModel, which="MBOprogress") <- optVal$MBOprogress
  return(finalModel)
}

########################################
# Reformulate EI without numerical check
# Check is flawed and gives error in some cases

EImod <- function (x, model, plugin = NULL, type = "UK", minimization = TRUE, 
                   envir = NULL) 
{
  if (is.null(plugin)) {
    if (minimization) {
      plugin <- min(model@y)
    }
    else {
      plugin <- -max(model@y)
    }
  }
  m <- plugin
  d <- length(x)
  if (d != model@d) {
    stop("x does not have the right size")
  }
  newdata.num <- as.numeric(x)
  newdata <- data.frame(t(newdata.num))
  colnames(newdata) = colnames(model@X)
  predx <- predict(object = model, newdata = newdata, type = type, 
                   checkNames = FALSE)
  kriging.mean <- predx$mean
  if (!minimization) {
    kriging.mean <- -kriging.mean
  }
  kriging.sd <- predx$sd
  xcr <- (m - kriging.mean)/kriging.sd
  #######################
  # Code in function EI()
  #   if (kriging.sd/sqrt(model@covariance@sd2) < 1e-06) {
  #     res <- 0
  #     xcr <- xcr.prob <- xcr.dens <- NULL
  #   }
  # else{
  #  xcr.prob <- pnorm(xcr)
  #  xcr.dens <- dnorm(xcr)
  #  res <- (m - kriging.mean) * xcr.prob + kriging.sd * xcr.dens
  # }
  
  xcr.prob <- pnorm(xcr)
  xcr.dens <- dnorm(xcr)
  res <- (m - kriging.mean) * xcr.prob + kriging.sd * xcr.dens
  
  # Numerical check of kriging.mean, kriging.sd inputs in EI formula
  CheckCondition <- is.infinite(res) | 
                    is.null(res) | 
                    is.nan(res) |
                    is.na(res)
  if(any(CheckCondition)) {
    res[CheckCondition] <- 0
  }
  
  if (!is.null(envir)) {
    assign("xcr", xcr, envir = envir)
    assign("xcr.prob", xcr.prob, envir = envir)
    assign("xcr.dens", xcr.dens, envir = envir)
    assign("kriging.sd", kriging.sd, envir = envir)
    assign("c", predx$c, envir = envir)
    assign("Tinv.c", predx$Tinv.c, envir = envir)
  }
  return(res)
}

#########################################
# MBO with shared hyperparameters + Level
#########################################

tuneMboSharedCvKDSN <- function (y, X, alphaShared=0, fineTuneIt=100, 
                                nStepMult=20, designMult=10, 
                                dimMax=round(sqrt(dim(X)[1])/2), addInfo=TRUE, # nCores=1,
                                maxRunsMult=1, repMult=1, tol_input=.Machine$double.eps^0.25, 
                                cvIndex=NULL, lossFunc=devStandard, EIopt="1Dmulti", GenSAmaxCall=100,
                                varSelectShared=FALSE, rdcRep=1, dropHiddenShared=FALSE,
                                standX=TRUE, standY=FALSE, timeAlloc="constant", varPreSelect=FALSE,
                                varPreSelpopSize=100, varPreSelMaxiter=100, maxLevels=10, 
                                useCV=TRUE, yTest=NULL, Xtest=NULL, EItype="EQI") {
  
#  if(nCores<1){stop("Please specify an integer number greater or equal to one as the number of threads!")}

  preSelectIndices <- NULL
  if(varPreSelect) {
    # Exclude columns prior variable selection with zero variance
    zeroSdVars <- which(colSds(X)==0)
    convertIndices <- setdiff(1:ncol(X), zeroSdVars)
    if(length(zeroSdVars)!=0) {
      X <- X[, -zeroSdVars, drop=FALSE]
    }

    # RDC variable selection based on all available covariates
    preSelectIndices <- rdcVarSelSubset(x=X, y=y, k=20, 
                                        s=1/6, f=sin, seedX=1:rdcRep, seedY=-c(1:rdcRep), 
                                        rdcRep=rdcRep, popSize=varPreSelpopSize, 
                                        maxiter=varPreSelMaxiter, addInfo=addInfo)

    # Select important variables
    X <- X[, preSelectIndices, drop=FALSE]
    
    # Convert indices to original scale
    preSelectIndices <- convertIndices[preSelectIndices]
    
    if(!useCV) {
      Xtest <- Xtest [, preSelectIndices, drop=FALSE]
    }
    if(addInfo) {cat("RDC variable pre-selection", "\n")}
  }
  
  rdcIndizes <- NULL
  if(varSelectShared) {
    rdcIndizes <- rdcVarOrder(x=X, y=y, cutoff=0, # nCores=nCores, 
                              seedX=1:rdcRep, seedY=-(1:rdcRep), rdcRep=rdcRep)
    if(addInfo) {cat("RDC marginal variable ordering", "\n")}
  }
  
  # Initialize parameters
  n <- dim(X)[1]
  
  # Initialize starting vector of hyperparameters
  levelsStart <- ifelse(round(maxLevels/2)!=0, round(maxLevels/2), 1)
  varExFreqStart <- 0.5
  MaxMatEntries <- .Machine$integer.max
  dimStart <- round ((dimMax+1) / 2)
  if((dimMax*2*n) > MaxMatEntries) {
    dimMax <- floor(MaxMatEntries/n/2)
    dimStart <- round (sqrt(dimMax*2)/2)
  }
  quantEuklid <- quantile(c(dist(robustStandard(X))^2), probs=c(0, 0.5, 1))
  # Correct for possible zero as lower bound of sigma
  if(any(quantEuklid==0)) {
    quantEuklid[quantEuklid==0] <- .Machine$double.eps^0.5
  }
  sigmaStart <- quantEuklid [2]
  lambdaRelStart <- 0.5
  dropHiddenSharedProbStart <- 0.5
  
  # Preconfigure tuning parameters and bounds
  noPar <- as.numeric(varSelectShared) + as.numeric(dropHiddenShared) + 4
  x_new <- vector("numeric", noPar)
  lower_bounds <- vector("numeric", noPar)
  upper_bounds <- vector("numeric", noPar)
  nameVec <- vector("character", noPar)
  if(varSelectShared & dropHiddenShared) {
    x_new <- c(levelsStart, varExFreqStart, dimStart, sigmaStart, dropHiddenSharedProbStart, lambdaRelStart)
    lower_bounds <- c(1, 0, 1, quantEuklid [1], 0.025, 0)
    upper_bounds <- c(maxLevels, 1, dimMax, quantEuklid [3], 0.975, 0.99)
    nameVec <- c("levels", "select", "dim", "sigma", "dropProb", "lambdaRel")
  }
  if(!varSelectShared & dropHiddenShared) {
    x_new <- c(levelsStart, dimStart, sigmaStart, dropHiddenSharedProbStart, lambdaRelStart)
    lower_bounds <- c(1, 1, quantEuklid [1], 0.025, 0)
    upper_bounds <- c(maxLevels, dimMax, quantEuklid [3], 0.975, 0.99)
    nameVec <- c("levels", "dim", "sigma", "dropProb", "lambdaRel")
  }
  if(varSelectShared & !dropHiddenShared) {
    x_new <- c(levelsStart, varExFreqStart, dimStart, sigmaStart, lambdaRelStart)
    lower_bounds <- c(1, 0, 1, quantEuklid [1], 0)
    upper_bounds <- c(maxLevels, 1, dimMax, quantEuklid [3], 0.99)
    nameVec <- c("levels", "select", "dim", "sigma", "lambdaRel")
  }
  if(!varSelectShared & !dropHiddenShared) {
    x_new <- c(levelsStart, dimStart, sigmaStart, lambdaRelStart)
    lower_bounds <- c(1, 1, quantEuklid [1], 0)
    upper_bounds <- c(maxLevels, dimMax, quantEuklid [3], 0.99)
    nameVec <- c("levels", "dim", "sigma", "lambdaRel")
  }

  if(useCV) {
  # Tune with cross validation
  optVal <- mboAll (loss_func=function (x) lossSharedCvKDSN (parOpt=x, y=y, X=X, cvIndex=cvIndex,  
                                                       seedW=rep(0, round(x[1])), lossFunc=lossFunc,
                                                       varSelectShared=varSelectShared, varRanking=rdcIndizes, alphaShared=alphaShared,
                                                       dropHiddenShared=dropHiddenShared, seedDrop=rep(0, round(x[1])), 
                                                       standX=standX, standY=standY, namesParOpt=nameVec), 
                    n_steps=(length(x_new)+3)*nStepMult, initDesign=(length(x_new)+3)*designMult, 
                    lower_bounds=lower_bounds, upper_bounds=upper_bounds, x_start=x_new,
                    # nCores=nCores, 
                    maxRunsMult=maxRunsMult, repMult=repMult, tol_input=tol_input, addInfo=addInfo,
                    EIopt=EIopt, GenSAmaxCall=GenSAmaxCall, timeAlloc=timeAlloc, EItype=EItype)
  }
  else{
    optVal <- mboAll (loss_func=function (x) lossSharedTestKDSN (parOpt=x, y=y, X=X, yTest=yTest, Xtest=Xtest,  
                                                               seedW=rep(0, round(x[1])), lossFunc=lossFunc,
                                                               varSelectShared=varSelectShared, varRanking=rdcIndizes, alphaShared=alphaShared,
                                                               dropHiddenShared=dropHiddenShared, seedDrop=rep(0, round(x[1])), 
                                                               standX=standX, standY=standY, namesParOpt=nameVec), 
                      n_steps=(length(x_new)+3)*nStepMult, initDesign=(length(x_new)+3)*designMult, 
                      lower_bounds=lower_bounds, upper_bounds=upper_bounds, x_start=x_new,
                      # nCores=nCores, 
                      maxRunsMult=maxRunsMult, repMult=repMult, tol_input=tol_input, addInfo=addInfo, 
                      EIopt=EIopt, GenSAmaxCall=GenSAmaxCall, timeAlloc=timeAlloc, EItype=EItype)
  }
  x_new <- optVal$par
  
  # Fine tune random fourier transformation weights
  # Reproduceability is ensured with seed generation
  if(fineTuneIt > 0) {
    localEnvir <- environment()
    fineTune <- vector("numeric", fineTuneIt)
    seedGuess <- matrix(sample.int(.Machine$integer.max, size=fineTuneIt * round(x_new[1])) * sample(c(-1, 1), 
                        size=fineTuneIt * round(x_new[1]), replace=TRUE), 
                        nrow=fineTuneIt, ncol=round(x_new[1]))
    seedDropGuess <- matrix(sample.int(.Machine$integer.max, size=fineTuneIt * round(x_new[1])) * sample(c(-1, 1), 
                            size=fineTuneIt * round(x_new[1]), replace=TRUE), 
                            nrow=fineTuneIt, ncol=round(x_new[1]))
#    if(nCores==1) {
      for(i in 1:fineTuneIt) {
        fineTune[i] <- lossSharedCvKDSN(parOpt=x_new, y=y, X=X, cvIndex=cvIndex, 
                                  seedW=seedGuess [i, ], lossFunc=lossFunc, 
                                  varSelectShared=varSelectShared, varRanking=rdcIndizes, 
                                  alphaShared=alphaShared,
                                  dropHiddenShared=dropHiddenShared,
                                  seedDrop=seedDropGuess[i, ], 
                                  standX=standX, standY=standY, namesParOpt=nameVec)[1]
        if(addInfo) {cat("tuneKDSN", "FineTune =", i, "\n")}
      }
#    }
    # else {
    #   cl <- makeCluster(nCores)
    #   clusterExport(cl = cl, varlist=c("lossSharedCvKDSN", "x_new", "y", "X", "seedGuess", "fineTuneIt"), 
    #                 envir = localEnvir)
    #   fineTune <- parSapply(cl=cl, X=1:fineTuneIt, 
    #                         FUN=function(i) lossSharedCvKDSN(parOpt=x_new, y=y, X=X, cvIndex=cvIndex,  
    #                                                    seedW=seedGuess [i, ], lossFunc=lossFunc,
    #                                                    varSelectShared=varSelectShared, varRanking=rdcIndizes,
    #                                                    alphaShared=alphaShared,
    #                                                    dropHiddenShared=dropHiddenShared,
    #                                                    seedDrop=seedDropGuess[i, ],
    #                                                    standX=standX, standY=standY, namesParOpt=nameVec)[1])
    #   stopCluster(cl=cl)
    #   if(addInfo) {cat("FineTuning", "\n")}
    # }
    minIndex <- which.min(fineTune)
  } 
  else {
    fineTune <- Inf
    minIndex <- 1
  }
  
  # Refit best model
  # Extract tuning parameters
  levels <- round(x_new[grep("levels", nameVec)])
  if(varSelectShared) {
    varExFreq <- rep(x_new[grep("select", nameVec)], levels)
  }
  else{
    varExFreq <- rep(NA, levels)
  }
  Dim <- rep(round(x_new[grep("dim", nameVec)]), levels)
  sigma <- rep(x_new[grep("sigma", nameVec)], levels)
  if(dropHiddenShared) {
    dropHiddenProb <- rep(x_new[grep("dropProb", nameVec)], levels)
  }
  else{
    dropHiddenProb <- rep(NA, levels)
  }
  lambdaRel <- rep(x_new[grep("lambdaRel", nameVec)], levels)
  # Fixed parameters
  alpha <- rep(alphaShared, levels)
  dropHidden <- rep(dropHiddenShared, levels)
  varSelect <- rep(varSelectShared, levels)
  
  # Output
  if(fineTune[minIndex] < optVal$value){
    finalModel <- fitKDSN(y = y, X = X, levels = levels, Dim = Dim, 
                          sigma = sigma, lambdaRel = lambdaRel, alpha = alpha, 
                          info = FALSE, seedW = seedGuess [minIndex, ], 
                          standX = standX, standY=standY,
                          varSelect=varSelect, varRanking=rdcIndizes,
                          varExFreq=varExFreq, 
                          dropHidden=dropHidden, 
                          dropHiddenProb=dropHiddenProb,
                          seedDrop=seedDropGuess[minIndex, ])
    attr(finalModel, which="Loss") <- fineTune[minIndex]
  }
  else{
    finalModel <- fitKDSN(y = y, X = X, levels = levels, Dim = Dim, 
                          sigma = sigma, lambdaRel = lambdaRel, alpha = alpha, 
                          info = FALSE, seedW = rep(0, levels), 
                          standX = standX, standY=standY,
                          varSelect=varSelect, varRanking=rdcIndizes,
                          varExFreq=varExFreq,
                          dropHidden=dropHidden, 
                          dropHiddenProb=dropHiddenProb,
                          seedDrop=rep(0, levels))
    attr(finalModel, which="Loss") <- optVal$value
  }
  # Include Loss score, loss function, cvIndex and preselected variables as attributes
  attr(finalModel, which="LossFunc") <- lossFunc
  attr(finalModel, which="cvIndex") <- cvIndex
  if(varPreSelect) {
    finalModel$Output$varPreSelect <- TRUE
    attr(finalModel, which="preSelectVars") <- preSelectIndices
  }
  # Include progress of MBO optimization
  dimnames(optVal$MBOprogress)[[2]] <- c("Loss", nameVec)
  attr(finalModel, which="MBOprogress") <- optVal$MBOprogress
  return(finalModel)
}

######################################
# Tune MBO SKDSN on training subsets 

tuneMboSharedSubsetKDSN <- function(noSubsets=2, noSubsetRep=1, subSeed=NULL, y, X, 
                                    alphaShared=1, 
                                    nStepMult=20, designMult=10,
                                    lossFunc=devStandard, GenSAmaxCall=100,
                                    varSelectShared=TRUE, dropHiddenShared=TRUE,
                                    standX=TRUE, standY=FALSE, timeAlloc="constant", 
                                    varPreSelect=TRUE, varPreSelpopSize=100, 
                                    varPreSelMaxiter=100, maxLevels=10, nCores=1,
                                    addInfo=1, saveOnDisk=FALSE,
                                    dirNameDisk=paste(tempdir(), "/ensembleModel", sep=""),
                                    useAllSub=TRUE, trainPartition=0.5, noPartion=1,
                                    EItype="EQI") {

  if(!saveOnDisk) {
    
    if(nCores==1) {
      
      # Tune MBO on all training subsets
      if(useAllSub) {
        tunedSKDSNsub <- vector("list", noSubsets*noSubsetRep)
        designGrid <- expand.grid(Subset=1:noSubsets, Rep=1:noSubsetRep)
      }
      else{
        tunedSKDSNsub <- vector("list", noPartion)
        designGrid <- expand.grid(Rep=1:noPartion)
      }
      noDesignRows <- nrow(designGrid)
      addInfoInner <- ifelse(addInfo==2, TRUE, FALSE)
      for(i in 1:noDesignRows) {
        
        if(useAllSub) {
          # Generate subsets
          set.seed(subSeed[ designGrid[i, "Rep"] ])
          subsetInd <- createFolds(y=y, k=noSubsets, returnTrain=FALSE)
          
          # MBO tuning
          tunedSKDSNsub[[i]] <- tuneMboSharedCvKDSN (y=y[ subsetInd[[ designGrid[i, "Subset"] ]], , drop=FALSE], 
                                                     X=X[ subsetInd[[ designGrid[i, "Subset"] ]], , drop=FALSE], 
                                                     alphaShared=alphaShared, fineTuneIt=0, 
                                                     nStepMult=nStepMult, designMult=designMult, 
                                                     dimMax=length(subsetInd[[ designGrid[i, "Subset"] ]]), 
                                                     addInfo=addInfoInner,
                                                     lossFunc=lossFunc, EIopt="GenSA", GenSAmaxCall=GenSAmaxCall,
                                                     varSelectShared=varSelectShared, rdcRep=1, 
                                                     dropHiddenShared=dropHiddenShared,
                                                     standX=standX, standY=standY, timeAlloc=timeAlloc, 
                                                     varPreSelect=varPreSelect,
                                                     varPreSelpopSize=varPreSelpopSize, 
                                                     varPreSelMaxiter=varPreSelMaxiter, 
                                                     maxLevels=maxLevels, 
                                                     useCV=FALSE, 
                                                     yTest=y[ -subsetInd[[ designGrid[i, "Subset"] ]] , , drop=FALSE], 
                                                     Xtest=X[-subsetInd[[ designGrid[i, "Subset"] ]], , drop=FALSE],
                                                     EItype=EItype)
        }
        else{
          set.seed(subSeed[ designGrid[i, "Rep"] ])
          subsetInd <- createDataPartition(y=y, times=noPartion, p=trainPartition, list=TRUE)
          
          # MBO tuning
          tunedSKDSNsub[[i]] <- tuneMboSharedCvKDSN (y=y[ subsetInd[[ designGrid[i, "Rep"] ]], , drop=FALSE], 
                                                     X=X[ subsetInd[[ designGrid[i, "Rep"] ]], , drop=FALSE], 
                                                     alphaShared=alphaShared, fineTuneIt=0, 
                                                     nStepMult=nStepMult, designMult=designMult, 
                                                     dimMax=length(subsetInd[[ designGrid[i, "Rep"] ]]), 
                                                     addInfo=addInfoInner,
                                                     lossFunc=lossFunc, EIopt="GenSA", GenSAmaxCall=GenSAmaxCall,
                                                     varSelectShared=varSelectShared, rdcRep=1, 
                                                     dropHiddenShared=dropHiddenShared,
                                                     standX=standX, standY=standY, timeAlloc=timeAlloc, 
                                                     varPreSelect=varPreSelect,
                                                     varPreSelpopSize=varPreSelpopSize, 
                                                     varPreSelMaxiter=varPreSelMaxiter, 
                                                     maxLevels=maxLevels, 
                                                     useCV=FALSE, 
                                                     yTest=y[ -subsetInd[[ designGrid[i, "Rep"] ]] , , drop=FALSE], 
                                                     Xtest=X[-subsetInd[[ designGrid[i, "Rep"] ]], , drop=FALSE],
                                                     EItype=EItype)
        }

        if(addInfo==1 | addInfo==2){
          cat("Subset", designGrid[i, "Subset"], 
              "Repetition", designGrid[i, "Rep"], 
              "Progress", round(i/noDesignRows, 4)*100, "%", "\n")
        }
      }
    }
    else{
      
      # Calculate grid
      if(useAllSub) {
        designGrid <- expand.grid(Subset=1:noSubsets, Rep=1:noSubsetRep)
      }
      else{
        designGrid <- expand.grid(Rep=1:noPartion)
      }
      noDesignRows <- nrow(designGrid)
      
      # Define temporary function
      tempFun <- function(i) {
        if(useAllSub) {
          # Generate subsets
          set.seed(subSeed[ designGrid[i, "Rep"] ])
          subsetInd <- createFolds(y=y, k=noSubsets, returnTrain=FALSE)
          
          # MBO
          tempFit <- tuneMboSharedCvKDSN (y=y[ subsetInd[[ designGrid[i, "Subset"] ]] , , drop=FALSE], 
                               X=X[ subsetInd[[ designGrid[i, "Subset"] ]], , drop=FALSE], 
                               alphaShared=alphaShared, fineTuneIt=0, 
                               nStepMult=nStepMult, designMult=designMult, 
                               dimMax=length(subsetInd[[ designGrid[i, "Subset"] ]]), 
                               addInfo=FALSE,
                               lossFunc=lossFunc, EIopt="GenSA", GenSAmaxCall=GenSAmaxCall,
                               varSelectShared=varSelectShared, rdcRep=1, 
                               dropHiddenShared=dropHiddenShared,
                               standX=standX, standY=standY, timeAlloc=timeAlloc, 
                               varPreSelect=varPreSelect,
                               varPreSelpopSize=varPreSelpopSize, 
                               varPreSelMaxiter=varPreSelMaxiter, 
                               maxLevels=maxLevels, 
                               useCV=FALSE, 
                               yTest=y[ -subsetInd[[ designGrid[i, "Subset"] ]] , , drop=FALSE], 
                               Xtest=X[-subsetInd[[ designGrid[i, "Subset"] ]], , drop=FALSE],
                               EItype=EItype)
          return(tempFit)
        }
        else{
          # Generate train partition
          set.seed(subSeed[ designGrid[i, "Rep"] ])
          subsetInd <- createDataPartition(y=y, times=noPartion, p=trainPartition, list=TRUE)
          
          # MBO
          tempFit <- tuneMboSharedCvKDSN (y=y[ subsetInd[[ designGrid[i, "Rep"] ]] , , drop=FALSE], 
                               X=X[ subsetInd[[ designGrid[i, "Rep"] ]], , drop=FALSE], 
                               alphaShared=alphaShared, fineTuneIt=0, 
                               nStepMult=nStepMult, designMult=designMult, 
                               dimMax=length(subsetInd[[ designGrid[i, "Rep"] ]]), 
                               addInfo=FALSE,
                               lossFunc=lossFunc, EIopt="GenSA", GenSAmaxCall=GenSAmaxCall,
                               varSelectShared=varSelectShared, rdcRep=1, 
                               dropHiddenShared=dropHiddenShared,
                               standX=standX, standY=standY, timeAlloc=timeAlloc, 
                               varPreSelect=varPreSelect,
                               varPreSelpopSize=varPreSelpopSize, 
                               varPreSelMaxiter=varPreSelMaxiter, 
                               maxLevels=maxLevels, 
                               useCV=FALSE, 
                               yTest=y[ -subsetInd[[ designGrid[i, "Rep"] ]] , , drop=FALSE], 
                               Xtest=X[-subsetInd[[ designGrid[i, "Rep"] ]], , drop=FALSE],
                               EItype=EItype)
          return(tempFit)
        }
        
      }
      
      # Cluster setup
      clustStart <- makeCluster(nCores)
      # Evaluate all variables in workspace (necessary because of lazy evaluation)
      eval(parse(text=ls()))
      # Export from current environment
      clusterExport(cl=clustStart, varlist=ls(), envir=environment())
      tunedSKDSNsub <- parLapply(cl=clustStart, X=1:noDesignRows, fun=tempFun)
    }
    
    # Output
    class(tunedSKDSNsub) <- "KDSNensemble"
    return(tunedSKDSNsub)
  }
  else{
    
    if(nCores==1) {
      
      # Tune MBO on all training subsets
      if(useAllSub) {
        tunedSKDSNsub <- vector("list", noSubsets*noSubsetRep)
        designGrid <- expand.grid(Subset=1:noSubsets, Rep=1:noSubsetRep)
      }
      else{
        tunedSKDSNsub <- vector("list", noPartion)
        designGrid <- expand.grid(Rep=1:noPartion)
      }
      noDesignRows <- nrow(designGrid)
      filePathVec <- vector("character", noDesignRows)
      addInfoInner <- ifelse(addInfo==2, TRUE, FALSE)
      for(i in 1:noDesignRows) {
        
        # Generate subsets
        if(useAllSub) {
          set.seed(subSeed[ designGrid[i, "Rep"] ])
          subsetInd <- createFolds(y=y, k=noSubsets, returnTrain=FALSE)
          
          # MBO tuning
          ensembleModel <- tuneMboSharedCvKDSN (y=y[ subsetInd[[ designGrid[i, "Subset"] ]] , , drop=FALSE], 
                                                X=X[ subsetInd[[ designGrid[i, "Subset"] ]], , drop=FALSE], 
                                                alphaShared=alphaShared, fineTuneIt=0, 
                                                nStepMult=nStepMult, designMult=designMult, 
                                                dimMax=length(subsetInd[[ designGrid[i, "Subset"] ]]), 
                                                addInfo=addInfoInner,
                                                lossFunc=lossFunc, EIopt="GenSA", GenSAmaxCall=GenSAmaxCall,
                                                varSelectShared=varSelectShared, rdcRep=1, 
                                                dropHiddenShared=dropHiddenShared,
                                                standX=standX, standY=standY, timeAlloc=timeAlloc, 
                                                varPreSelect=varPreSelect,
                                                varPreSelpopSize=varPreSelpopSize, 
                                                varPreSelMaxiter=varPreSelMaxiter, 
                                                maxLevels=maxLevels, 
                                                useCV=FALSE, 
                                                yTest=y[ -subsetInd[[ designGrid[i, "Subset"] ]] , , drop=FALSE], 
                                                Xtest=X[-subsetInd[[ designGrid[i, "Subset"] ]], , drop=FALSE],
                                                EItype=EItype)
        }
        else{
          set.seed(subSeed[ designGrid[i, "Rep"] ])
          subsetInd <- createDataPartition(y=y, times=noPartion, p=trainPartition, list=TRUE)
          
          # MBO tuning
          ensembleModel <- tuneMboSharedCvKDSN (y=y[ subsetInd[[ designGrid[i, "Rep"] ]] , , drop=FALSE], 
                                                X=X[ subsetInd[[ designGrid[i, "Rep"] ]], , drop=FALSE], 
                                                alphaShared=alphaShared, fineTuneIt=0, 
                                                nStepMult=nStepMult, designMult=designMult, 
                                                dimMax=length(subsetInd[[ designGrid[i, "Rep"] ]]), 
                                                addInfo=addInfoInner,
                                                lossFunc=lossFunc, EIopt="GenSA", GenSAmaxCall=GenSAmaxCall,
                                                varSelectShared=varSelectShared, rdcRep=1, 
                                                dropHiddenShared=dropHiddenShared,
                                                standX=standX, standY=standY, timeAlloc=timeAlloc, 
                                                varPreSelect=varPreSelect,
                                                varPreSelpopSize=varPreSelpopSize, 
                                                varPreSelMaxiter=varPreSelMaxiter, 
                                                maxLevels=maxLevels, 
                                                useCV=FALSE, 
                                                yTest=y[ -subsetInd[[ designGrid[i, "Rep"] ]] , , drop=FALSE], 
                                                Xtest=X[-subsetInd[[ designGrid[i, "Rep"] ]], , drop=FALSE],
                                                EItype=EItype)
        }

        # Save model on disk and remove it from workspace
        filePathVec[i] <- paste(dirNameDisk, i, sep="")
        save(ensembleModel, file=filePathVec[i], compress="xz")
        rm(ensembleModel)
        
        if(addInfo==1 | addInfo==2){
          cat("Subset", designGrid[i, "Subset"], 
              "Repetition", designGrid[i, "Rep"], 
              "Progress", round(i/noDesignRows, 4)*100, "%", "\n")
        }
      }
    }
    else{
      
      # Calculate grid
      if(useAllSub) {
        designGrid <- expand.grid(Subset=1:noSubsets, Rep=1:noSubsetRep)
      }
      else{
        designGrid <- expand.grid(Rep=1:noPartion)
      }
      noDesignRows <- nrow(designGrid)
      
      # Define temporary function
      tempFun <- function(i) {
        # Generate subsets
        if(useAllSub) {
          set.seed(subSeed[ designGrid[i, "Rep"] ])
          subsetInd <- createFolds(y=y, k=noSubsets, returnTrain=FALSE)
          
          ensembleModel <- tuneMboSharedCvKDSN (y=y[ subsetInd[[ designGrid[i, "Subset"] ]] , , drop=FALSE], 
                                                X=X[subsetInd[[ designGrid[i, "Subset"] ]], , drop=FALSE], 
                                                alphaShared=alphaShared, fineTuneIt=0, 
                                                nStepMult=nStepMult, designMult=designMult, 
                                                dimMax=length(subsetInd[[ designGrid[i, "Subset"] ]]), addInfo=FALSE,
                                                lossFunc=lossFunc, EIopt="GenSA", GenSAmaxCall=GenSAmaxCall,
                                                varSelectShared=varSelectShared, rdcRep=1, 
                                                dropHiddenShared=dropHiddenShared,
                                                standX=standX, standY=standY, timeAlloc=timeAlloc, 
                                                varPreSelect=varPreSelect,
                                                varPreSelpopSize=varPreSelpopSize, 
                                                varPreSelMaxiter=varPreSelMaxiter, 
                                                maxLevels=maxLevels, 
                                                useCV=FALSE, 
                                                yTest=y[ -subsetInd[[ designGrid[i, "Subset"] ]] , , drop=FALSE], 
                                                Xtest=X[-subsetInd[[ designGrid[i, "Subset"] ]], , drop=FALSE],
                                                EItype=EItype)
        }
        else{
          set.seed(subSeed[ designGrid[i, "Rep"] ])
          subsetInd <- createDataPartition(y=y, times=noPartion, p=trainPartition, list=TRUE)
          
          ensembleModel <- tuneMboSharedCvKDSN (y=y[ subsetInd[[ designGrid[i, "Rep"] ]] , , drop=FALSE], 
                                                X=X[subsetInd[[ designGrid[i, "Rep"] ]], , drop=FALSE], 
                                                alphaShared=alphaShared, fineTuneIt=0, 
                                                nStepMult=nStepMult, designMult=designMult, 
                                                dimMax=length(subsetInd[[ designGrid[i, "Rep"] ]]), addInfo=FALSE,
                                                lossFunc=lossFunc, EIopt="GenSA", GenSAmaxCall=GenSAmaxCall,
                                                varSelectShared=varSelectShared, rdcRep=1, 
                                                dropHiddenShared=dropHiddenShared,
                                                standX=standX, standY=standY, timeAlloc=timeAlloc, 
                                                varPreSelect=varPreSelect,
                                                varPreSelpopSize=varPreSelpopSize, 
                                                varPreSelMaxiter=varPreSelMaxiter, 
                                                maxLevels=maxLevels, 
                                                useCV=FALSE, 
                                                yTest=y[ -subsetInd[[ designGrid[i, "Rep"] ]] , , drop=FALSE], 
                                                Xtest=X[-subsetInd[[ designGrid[i, "Rep"] ]], , drop=FALSE],
                                                EItype=EItype)
        }
        
        # Save model on disk and remove it from workspace
        filePath <- paste(dirNameDisk, i, sep="")
        save(ensembleModel, file=filePath, compress="xz")
        rm(ensembleModel)
        return(filePath)
      }
      
      # Cluster setup
      clustStart <- makeCluster(nCores)
      # Evaluate all variables in workspace (necessary because lazy evaluation)
      eval(parse(text=ls()))
      # Export from current environment
      clusterExport(cl=clustStart, varlist=ls(), envir=environment())
      filePathVec <- parSapply(cl=clustStart, X=1:noDesignRows, FUN=tempFun)
    }
    
    # Output
    class(filePathVec) <- "KDSNensembleDisk"
    return(filePathVec)
  }
}

##############################################
# Tune with GCV approximation without test set
##############################################

tuneMboLevelGcvKDSN <- function (y, X, levels=1, alpha=rep(0, levels), fineTuneIt=100, 
                                nStepMult=20, designMult=10, 
                                dimMax=round(sqrt(dim(X)[1])/2), addInfo=TRUE,
                                maxRunsMult=1, repMult=1, tol_input=.Machine$double.eps^0.25, 
                                EIopt="1Dmulti", GenSAmaxCall=100,
                                varSelect=rep(FALSE, levels), rdcRep=1, dropHidden=rep(FALSE, levels),
                                standX=TRUE, standY=FALSE, timeAlloc="constant", varPreSelect=FALSE,
                                varPreSelpopSize=100, varPreSelMaxiter=100, EItype="EQI") {

  preSelectIndices <- NULL
  if(varPreSelect) {
    # Exclude columns prior variable selection with zero variance
    zeroSdVars <- which(colSds(X)==0)
    convertIndices <- setdiff(1:ncol(X), zeroSdVars)
    if(length(zeroSdVars)!=0) {
      X <- X[, -zeroSdVars, drop=FALSE]
    }
    
    # RDC variable selection based on all available covariates
    preSelectIndices <- rdcVarSelSubset(x=X, y=y, k=20, 
                                        s=1/6, f=sin, seedX=1:rdcRep, seedY=-c(1:rdcRep), 
                                        rdcRep=rdcRep, popSize=varPreSelpopSize, 
                                        maxiter=varPreSelMaxiter)
    
    # Select important variables
    X <- X[, preSelectIndices, drop=FALSE]
    
    # Convert indices to original scale for future predictions
    preSelectIndices <- convertIndices[preSelectIndices]
    
    if(addInfo) {cat("RDC variable pre-selection", "\n")}
  }
  
  rdcIndizes <- NULL
  if(any(varSelect)) {
    rdcIndizes <- rdcVarOrder(x=X, y=y, cutoff=0, # nCores=nCores, 
                              seedX=1:rdcRep, seedY=-(1:rdcRep), rdcRep=rdcRep)
    if(addInfo) {cat("RDC marginal variable ordering", "\n")}
  }
  
  # Initialize parameters
  n <- dim(X)[1]
  
  # Initialize starting vector of hyperparameters
  varExFreqStart <- 0.5
  MaxMatEntries <- .Machine$integer.max
  dimStart <- round ((dimMax+1) / 2)
  if((dimMax*2*n) > MaxMatEntries) {
    dimMax <- floor(MaxMatEntries/n/2)
    dimStart <- round (sqrt(dimMax*2)/2)
  }
  quantEuklid <- quantile(c(dist(robustStandard(X))^2), probs=c(0, 0.5, 1))
  # Correct for possible zero as lower bound of sigma
  if(any(quantEuklid==0)) {
    quantEuklid[quantEuklid==0] <- .Machine$double.eps^0.5
  }
  sigmaStart <- quantEuklid [2]
  lambdaRelStart <- 0.5
  dropHiddenProbStart <- 0.5

  # Consider four different cases:
  # 1: Variable selection and dropout
  # 2: Dropout
  # 3: Variable selection
  # 4: Neither variable selection nor dropout
  
  tempInd <- as.numeric(varSelect) + as.numeric(dropHidden) + 3
  noPar <- sum(tempInd)
  cumIndex <- c(0, cumsum(tempInd))
  Indices <- lapply(2:length(cumIndex), function(x) (1+cumIndex[x-1]):(cumIndex[x]))
  x_new <- vector("numeric", noPar)
  lower_bounds <- vector("numeric", noPar)
  upper_bounds <- vector("numeric", noPar)
  nameVec <- vector("character", noPar)
  for(l in 1:levels) {
    if(varSelect[l] & dropHidden[l]){
      x_new[Indices[[l]] ] <- c(varExFreqStart, dimStart, sigmaStart, dropHiddenProbStart, lambdaRelStart)
      lower_bounds[Indices[[l]] ] <- c(0, 1, quantEuklid [1], 0.025, 0)
      upper_bounds[Indices[[l]] ] <- c(1, dimMax, quantEuklid [3], 0.975, 0.99)
      nameVec[Indices[[l]] ] <- c("select", "dim", "sigma", "dropProb", "lambdaRel")
    }
    if(!varSelect[l] & dropHidden[l]){
      x_new[Indices[[l]] ] <- c(dimStart, sigmaStart, dropHiddenProbStart, lambdaRelStart)
      lower_bounds[Indices[[l]] ] <- c(1, quantEuklid [1], 0.025, 0)
      upper_bounds[Indices[[l]] ] <- c(dimMax, quantEuklid [3], 0.975, 0.99)
      nameVec[Indices[[l]] ] <- c("dim", "sigma", "dropProb", "lambdaRel")
    }
    if(varSelect[l] & !dropHidden[l]){
      x_new[Indices[[l]] ] <- c(varExFreqStart, dimStart, sigmaStart, lambdaRelStart)
      lower_bounds[Indices[[l]] ] <- c(0, 1, quantEuklid [1], 0)
      upper_bounds[Indices[[l]] ] <- c(1, dimMax, quantEuklid [3], 0.99)
      nameVec[Indices[[l]] ] <- c("select", "dim", "sigma", "lambdaRel")
    }
    if(!varSelect[l] & !dropHidden[l]){
      x_new[Indices[[l]] ] <- c(dimStart, sigmaStart, lambdaRelStart)
      lower_bounds[Indices[[l]] ] <- c(1, quantEuklid [1], 0)
      upper_bounds[Indices[[l]] ] <- c(dimMax, quantEuklid [3], 0.99)
      nameVec[Indices[[l]] ] <- c("dim", "sigma", "lambdaRel")
    }
  }

  # Tune
  optVal <- mboAll (loss_func= function (x) lossApprox (parOpt=x, y=y, X=X, levels=levels, seedW=seq(0, (levels-1), 1),
                                                        varSelect=varSelect, varRanking=rdcIndizes, 
                                                        alpha=alpha, dropHidden=dropHidden,
                                                        seedDrop=seq(0, -(levels-1), -1), standX=standX, standY=standY,
                                                        namesParOpt=nameVec,
                                                        gammaPar=1), 
                    n_steps=(length(x_new)+3)*nStepMult, initDesign=(length(x_new)+3)*designMult, 
                    lower_bounds=lower_bounds, upper_bounds=upper_bounds, x_start=x_new,
                    # nCores=nCores, 
                    addInfo=addInfo, 
                    maxRunsMult=maxRunsMult, repMult=repMult, tol_input=tol_input, 
                    EIopt=EIopt, GenSAmaxCall=GenSAmaxCall, timeAlloc=timeAlloc,
                    EItype=EItype)
  x_new <- optVal$par
  
  # Fine tune random fourier transformation weights
  # Reproduceability is ensured with seed generation
  if(fineTuneIt > 0) {
    fineTune <- vector("numeric", fineTuneIt)
    seedGuess <- matrix(sample.int(.Machine$integer.max, size=fineTuneIt * levels) * sample(c(-1, 1), size=fineTuneIt * levels, replace=TRUE), 
                        nrow=fineTuneIt, ncol=levels)
    seedDropGuess <- matrix(sample.int(.Machine$integer.max, size=fineTuneIt * levels) * sample(c(-1, 1), size=fineTuneIt * levels, replace=TRUE), 
                            nrow=fineTuneIt, ncol=levels)
    localEnvir <- environment()
    for(i in 1:fineTuneIt) {
      fineTune[i] <-  lossApprox(parOpt=x_new, y=y, X=X, levels=levels, seedW=seedGuess [i, ],
                                 varSelect=varSelect, varRanking=rdcIndizes, 
                                 alpha=alpha, dropHidden=dropHidden,
                                 seedDrop=seedDropGuess[i, ], standX=standX, standY=standY,
                                 namesParOpt=nameVec, gammaPar=1)[1]
      
      if(addInfo) {cat("tuneKDSN", "FineTune =", i, "\n")}
    }
    minIndex <- which.min(fineTune)
  }
  else {
    fineTune <- Inf
    minIndex <- 1
  }
  
  # Output
  # Refit best model
  varExFreq <- x_new[grep("select", nameVec)]
  Dim <- round(x_new[grep("dim", nameVec)])
  sigma <- x_new[grep("sigma", nameVec)]
  dropHiddenProb <- x_new[grep("dropProb", nameVec)]
  lambdaRel <- x_new[grep("lambdaRel", nameVec)]
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
  
  if(fineTune[minIndex] < optVal$value){
    finalModel <- fitKDSN(y = y, X = X, levels = levels, Dim = Dim, 
                          sigma = sigma, lambdaRel = lambdaRel, alpha = alpha, 
                          info = FALSE, seedW = seedGuess [minIndex, ], 
                          standX = standX, standY=standY,
                          varSelect=varSelect, varRanking=rdcIndizes,
                          varExFreq=varExFreq, 
                          dropHidden=dropHidden, 
                          dropHiddenProb=dropHiddenProb,
                          seedDrop=seedDropGuess[minIndex, ])
    attr(finalModel, which="Loss") <- fineTune[minIndex]
  }
  else{
    finalModel <- fitKDSN(y = y, X = X, levels = levels, Dim = Dim, 
                          sigma = sigma, lambdaRel = lambdaRel, alpha = alpha, 
                          info = FALSE, seedW = seq(0, (levels-1), 1), 
                          standX = standX, standY=standY,
                          varSelect=varSelect, varRanking=rdcIndizes,
                          varExFreq=varExFreq,
                          dropHidden=dropHidden, 
                          dropHiddenProb=dropHiddenProb,
                          seedDrop=seq(0, -(levels-1), -1))
    attr(finalModel, which="Loss") <- optVal$value
  }
  # Include Loss score, loss function, cvIndex and preselected variables as attributes
  attr(finalModel, which="LossFunc") <- "GCV"
  if(varPreSelect) {
    finalModel$Output$varPreSelect <- TRUE
    attr(finalModel, which="preSelectVars") <- preSelectIndices
  }
  # Include progress of MBO optimization
  dimnames(optVal$MBOprogress)[[2]] <- c("Loss", nameVec)
  attr(finalModel, which="MBOprogress") <- optVal$MBOprogress
  return(finalModel)
}
