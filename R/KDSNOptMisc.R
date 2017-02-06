################################
# Finetuning of given KDSN model

# With cross validation
fineTuneCvKDSN <- function(estKDSN, y, X, fineTuneIt=100, info=TRUE, 
                         seedInitW=NULL, seedInitD=NULL, ncpus=1,
                         cvIndex, lossFunc=devStandard) {
  # Input adjustments
  if(ncpus<1 | floor(ncpus)!=ncpus) {stop("Please supply a positive number of available cpus greater 
                                          or equal to 1 (integer scalar)")}
  fineTune <- vector("numeric", fineTuneIt)
  levels <- estKDSN$Input$levels
  x_new <- c(matrix(c(estKDSN$Input$varExFreq, estKDSN$Input$Dim, estKDSN$Input$sigma, 
                      estKDSN$Input$dropHiddenProb, estKDSN$Input$lambdaRel), 
                    byrow=TRUE, nrow=5))
  names_x_new <- rep(c("select", "dim", "sigma", "dropProb", "lambdaRel"), levels)
  if(sum(which(is.na(x_new)))!=0) {
    names_x_new <- names_x_new[-which(is.na(x_new))]
  }
  # Exclude values, which will not be optimized
  x_new <- na.omit(x_new)
  seedOriginalW <- estKDSN$Input$seedW
  seedOriginalD <- estKDSN$Input$seedDrop
  if(is.null(seedOriginalW)) {
#    seedOriginalW <- seq(0, (levels-1), 1)
    seedOriginalW <- rep(0, (levels-1))
  }
  if(is.null(seedOriginalD)) {
#    seedOriginalD <- seq(0, -(levels-1), -1)
    seedOriginalD <- rep(0, (levels-1))
  }

  # Reproduceability is ensured with seed generation
  set.seed(seedInitW)
  seedGuessW <- matrix(sample.int(.Machine$integer.max, size=fineTuneIt * levels) * sample(c(-1, 1), size=fineTuneIt * levels, replace=TRUE), nrow=fineTuneIt, ncol=levels)
  set.seed(seedInitD)
  seedGuessD <- matrix(sample.int(.Machine$integer.max, size=fineTuneIt * levels) * sample(c(-1, 1), size=fineTuneIt * levels, replace=TRUE), nrow=fineTuneIt, ncol=levels)
  envir1 <- environment()

  # Apply variabel pre-selection
  if(any(names(attributes(estKDSN))=="preSelectVars")){
    X <- X[, attr(estKDSN, "preSelectVars"), drop=FALSE]
  }
  
  if(ncpus==1) {
    # Fine tune random fourier transformation weights
    for(i in 1:fineTuneIt) {
      fineTune [i] <- lossCvKDSN (parOpt=x_new, y=y, X=X, levels=levels, seedW=seedGuessW [i, ],
                                cvIndex=cvIndex, lossFunc=lossFunc, alpha=estKDSN$Input$alpha, 
                                varSelect=estKDSN$Input$varSelect, varRanking=estKDSN$Input$varRanking,
                                dropHidden=estKDSN$Input$dropHidden, seedDrop=seedGuessD [i, ], 
                                namesParOpt=names_x_new)
      if(info) {cat("tuneKDSN", "FineTune =", i, "\n")}
    }
  }
  else {
    tempFunc <- function(i) {
      fineTuneResult <- lossCvKDSN (parOpt=x_new, y=y, X=X, levels=levels, seedW=seedGuessW [i, ],
                                    cvIndex=cvIndex, lossFunc=lossFunc, alpha=estKDSN$Input$alpha, 
                                    varSelect=estKDSN$Input$varSelect, varRanking=estKDSN$Input$varRanking,
                                    dropHidden=estKDSN$Input$dropHidden, seedDrop=seedGuessD [i, ], 
                                    namesParOpt=names_x_new)
      return(fineTuneResult)
    }
    
    # Parallel execution with parallel package
    cl <- makeCluster(ncpus)
    clusterExport(cl=cl, varlist=c("seedGuessW", "seedGuessD", "x_new", "y", "X", 
                                   "cvIndex", "lossFunc", "levels", "estKDSN", "tempFunc", 
                                   "names_x_new"), 
                  envir=envir1)
    clusterEvalQ(cl=cl, library(kernDeepStackNet))
    fineTune <- parSapply(cl=cl, X=1:fineTuneIt, FUN=tempFunc)
    stopCluster(cl)
  }
  minIndex <- which.min(fineTune)
  
  # Evaluate loss of given model
  givenLoss <- c(lossCvKDSN (parOpt=x_new, y=y, X=X, levels=levels, seedW=seedOriginalW,
                             cvIndex=cvIndex, lossFunc=lossFunc, alpha=estKDSN$Input$alpha, 
                             varSelect=estKDSN$Input$varSelect, varRanking=estKDSN$Input$varRanking,
                             dropHidden=estKDSN$Input$dropHidden, seedDrop=seedOriginalD, 
                             namesParOpt=names_x_new) )

  # Refit best model and output if it is better than the given model
  lenx_new <- length(x_new)
  stopifnot((lenx_new %% levels) == 0)
  stopifnot(length(seedGuessW[minIndex, ]) == levels)
  stopifnot(length(seedGuessD[minIndex, ]) == levels)
#   varExFreq1 <- x_new[seq(1, lenx_new, 5)]
#   Dim1 <- round(x_new[seq(2, lenx_new, 5)])
#   sigma1 <- x_new[seq(3, lenx_new, 5)]
#   dropHiddenProb1 <- x_new[seq(4, lenx_new, 5)]
#   lambdaRel1 <- x_new[seq(5, lenx_new, 5)]
  if(fineTune[minIndex] < givenLoss){
    finalModel <- fitKDSN(y = y, X = X, levels = levels, Dim = estKDSN$Input$Dim, 
                          sigma = estKDSN$Input$sigma, lambdaRel = estKDSN$Input$lambdaRel, 
                          alpha = estKDSN$Input$alpha,
                          varSelect=estKDSN$Input$varSelect, 
                          varExFreq=estKDSN$Input$varExFreq,
                          varRanking=estKDSN$Input$varRanking, 
                          dropHidden=estKDSN$Input$dropHidden,
                          dropHiddenProb=estKDSN$Input$dropHiddenProb, 
                          info = FALSE, 
                          seedW = seedGuessW [minIndex, ], 
                          standX = estKDSN$Input$standX,
                          standY = estKDSN$Input$standY,
                          seedDrop = seedGuessD [minIndex, ])
    # Include loss score as attribute
    attr(finalModel, which = "Loss") <- fineTune[minIndex]
    attr(finalModel, which = "LossFunc") <- lossFunc
    attr(finalModel, which = "cvIndex") <- cvIndex
    
    # VarPreSelect
    if (any(names(attributes(estKDSN))=="preSelectVars")) {
      attr(finalModel, which = "preSelectVars") <- attr(estKDSN, "preSelectVars")
    }
    return(finalModel)
  }
  else{
    return(estKDSN)
  }
}

##############################
# Ensemble based on KDSN model

fitEnsembleKDSN <- function(estKDSN, y, X, ensembleSize=100, 
                            bagging=rep(FALSE, ensembleSize), seedBag=NULL, 
                            randSubsets=rep(FALSE, ensembleSize), 
                            seedRandSubSize=NULL, seedRandSubDraw=NULL,
                            seedW1=NULL, seedW2=NULL, 
                            seedDrop1=NULL, seedDrop2=NULL, 
                            info=TRUE, nCores=1, saveOnDisk=FALSE, 
                            dirNameDisk=paste(tempdir(), "/ensembleModel", sep="")) {
  # Check input arguments
  stopifnot(class(estKDSN)=="KDSN")
  
  # Define fixed inputs
  levelsInput <- estKDSN$Input$levels
  DimInput <- estKDSN$Input$Dim
  sigmaInput <- estKDSN$Input$sigma
  lambdaRelInput <- estKDSN$Input$lambdaRel
  alphaInput <- estKDSN$Input$alpha
  standXinput <- estKDSN$Input$standX
  standYinput <- estKDSN$Input$standY
  varSelectInput <- estKDSN$Input$varSelect
  varRankingInput <- estKDSN$Input$varRanking
  varExFreqInput <- estKDSN$Input$varExFreq
  dropHiddenInput <- estKDSN$Input$dropHidden
  dropHiddenProbInput <- estKDSN$Input$dropHiddenProb
  varPreSelectInput <- any(names(attributes(estKDSN))=="preSelectVars")

  if(!saveOnDisk) {
  
  if(nCores==1) {
  
    # Variable pre selection
    if(varPreSelectInput) {
      Xinput <- X[, attr(estKDSN, "preSelectVars"), drop=FALSE]
    }
    else{
      Xinput <- X
    }
    
    # Fit all models
    ensembleModels <- vector("list", ensembleSize)
    for(i in 1:ensembleSize) {
      
      # Bagging
      if(bagging[i]) {
        # Draw random bootstrap samples
        set.seed(seedBag[i])
        bootSamples <- sample(1:length(y), size=levelsInput*length(y), replace=TRUE)
        # Convert to list
        bootSamples <- split(bootSamples, rep(1:levelsInput, each=length(y)))
      }
      else{
        bootSamples <- NULL
      }

      # Draw random Fourier transformation weights
      set.seed(seedW1[i])
      posNumbers <- sample.int(.Machine$integer.max, size=levelsInput)
      set.seed(seedW2[i])
      signNumbers <- sample(c(-1, 1), size=levelsInput, replace=TRUE)
      seedrandFourInput <- signNumbers * posNumbers
      
      # Draw random dropouts
      set.seed(seedDrop1[i])
      posNumbDrop <- sample.int(.Machine$integer.max, size=levelsInput)
      set.seed(seedDrop2[i])
      signNumbersDrop <- sample(c(-1, 1), size=levelsInput, replace=TRUE)
      seedDropInput <- signNumbersDrop * posNumbDrop
      
      # Random subsets
      if(randSubsets[i]) {
        
        # Draw size of random subset
        set.seed(seedRandSubSize[i])
        randSize <- sample(1:ncol(Xinput), size=levelsInput, replace=TRUE)
        
        # Draw random covariables given size as list
        randSubInd <- lapply(1:levelsInput, function(x) {
          seedInput <- ifelse(is.null(seedRandSubDraw[i]), NULL, seedRandSubDraw[i] + x)
          set.seed(seedInput)
          sample(1:ncol(Xinput), size=randSize[x], replace=FALSE)
        }
        )
      }
      else{
        randSubInd <- NULL
      }

      # Fit KDSN
      ensembleModels[[i]] <- fitKDSN(y = y, X = Xinput, levels = levelsInput, 
                                     Dim = DimInput, sigma = sigmaInput, 
                                     lambdaRel = lambdaRelInput, alpha = alphaInput, 
                                     info = FALSE, seedW = seedrandFourInput, 
                                     standX = standXinput, standY=standYinput,
                                     varSelect=varSelectInput, varRanking=varRankingInput,
                                     varExFreq=varExFreqInput, dropHidden=dropHiddenInput, 
                                     dropHiddenProb=dropHiddenProbInput,
                                     seedDrop=seedDropInput,
                                     baggingInd=bootSamples, randSubsetInd=randSubInd)
      
      # Include variable pre selection information, if applicable
      if(varPreSelectInput) {
        attr(ensembleModels[[i]], "preSelectVars") <- attr(estKDSN, "preSelectVars")
      }
      if(info) {
        cat("Ensemble", round(i/ensembleSize, 4)*100, "%", "\n")
      }
      
    }
    
    # Save additional inputs
    class(ensembleModels) <- "KDSNensemble"
    # Save bagging parameters
    attr(ensembleModels, "bagging") <- bagging
    attr(ensembleModels, "seedBag") <- seedBag
    # Save random subsets parameters
    attr(ensembleModels, "randSubsets") <- randSubsets
    attr(ensembleModels, "seedRandSubSize") <- seedRandSubSize
    attr(ensembleModels, "seedRandSubDraw") <- seedRandSubDraw
    # Save random fourier transformation seeds
    attr(ensembleModels, "seedW1") <- seedW1
    attr(ensembleModels, "seedW2") <- seedW2
    # Save random dropout seeds
    attr(ensembleModels, "seedDrop1") <- seedDrop1
    attr(ensembleModels, "seedDrop2") <- seedDrop2
    return(ensembleModels)
  }
  
  if(nCores>1) {
    localEnvir <- environment()
    tempFun <- function(i) {
      
      # Variable pre selection
      if(varPreSelectInput) {
        Xinput <- X[, attr(estKDSN, "preSelectVars")]
      }
      else{
        Xinput <- X
      }
      
      # Bagging
      if(bagging[i]) {
        # Draw random bootstrap samples
        set.seed(seedBag[i])
        bootSamples <- sample(1:length(y), size=levelsInput*length(y), replace=TRUE)
        # Convert to list
        bootSamples <- split(bootSamples, rep(1:levelsInput, each=length(y)))
      }
      else{
        bootSamples <- NULL
      }
      
      # Draw random Fourier transformation weights
      set.seed(seedW1[i])
      posNumbers <- sample.int(.Machine$integer.max, size=levelsInput)
      set.seed(seedW2[i])
      signNumbers <- sample(c(-1, 1), size=levelsInput, replace=TRUE)
      seedrandFourInput <- signNumbers * posNumbers
      
      # Draw random dropouts
      set.seed(seedDrop1[i])
      posNumbDrop <- sample.int(.Machine$integer.max, size=levelsInput)
      set.seed(seedDrop2[i])
      signNumbersDrop <- sample(c(-1, 1), size=levelsInput, replace=TRUE)
      seedDropInput <- signNumbersDrop * posNumbDrop
      
      # Random subsets
      if(randSubsets[i]) {
        
        # Draw size of random subset
        set.seed(seedRandSubSize[i])
        randSize <- sample(1:ncol(Xinput), size=levelsInput, replace=TRUE)
        
        # Draw random covariables given size as list
        randSubInd <- lapply(1:levelsInput, function(x) {
          seedInput <- ifelse(is.null(seedRandSubDraw[i]), NULL, seedRandSubDraw[i] + x)
          set.seed(seedInput)
          sample(1:ncol(Xinput), size=randSize[x], replace=FALSE)
        })
      }
      else{
        randSubInd <- NULL
      }
      
      # Fit KDSN
      ensembleModel <- fitKDSN(y = y, X = Xinput, levels = levelsInput, 
                               Dim = DimInput, sigma = sigmaInput, 
                               lambdaRel = lambdaRelInput, alpha = alphaInput, 
                               info = FALSE, seedW = seedrandFourInput, 
                               standX = standXinput, standY=standYinput,
                               varSelect=varSelectInput, varRanking=varRankingInput,
                               varExFreq=varExFreqInput, dropHidden=dropHiddenInput, 
                               dropHiddenProb=dropHiddenProbInput,
                               seedDrop=seedDropInput,
                               baggingInd=bootSamples, randSubsetInd=randSubInd)
      
      # Include variable pre selection information, if applicable
      if(varPreSelectInput) {
        attr(ensembleModel, "preSelectVars") <- attr(estKDSN, "preSelectVars")
      }
      return(ensembleModel)
    }
    parClust <- makeCluster(nCores)
    clusterExport(cl = parClust, varlist=c("levelsInput", "DimInput", "sigmaInput","lambdaRelInput", 
                                     "alphaInput", "standXinput", "standYinput", "varSelectInput",
                                     "varRankingInput", "varExFreqInput", "dropHiddenInput",
                                     "dropHiddenProbInput", "y", "X", "info",
                                     "bagging", "seedBag",
                                     "randSubsets", "seedRandSubSize", "seedRandSubDraw", 
                                     "seedW1", "seedW2", "varPreSelectInput",
                                     "seedDrop1", "seedDrop2"), envir = localEnvir)
    ensembleModels <- parLapply(cl=parClust, X=1:ensembleSize, fun=tempFun)
    stopCluster(parClust)
    
    # Save additional inputs
    class(ensembleModels) <- "KDSNensemble"
    # Save bagging parameters
    attr(ensembleModels, "bagging") <- bagging
    attr(ensembleModels, "seedBag") <- seedBag
    # Save random subsets parameters
    attr(ensembleModels, "randSubsets") <- randSubsets
    attr(ensembleModels, "seedRandSubSize") <- seedRandSubSize
    attr(ensembleModels, "seedRandSubDraw") <- seedRandSubDraw
    # Save random fourier transformation seeds
    attr(ensembleModels, "seedW1") <- seedW1
    attr(ensembleModels, "seedW2") <- seedW2
    # Save random dropout seeds
    attr(ensembleModels, "seedDrop1") <- seedDrop1
    attr(ensembleModels, "seedDrop2") <- seedDrop2
    return(ensembleModels)
  }
  }
  else{
    
    # Variable pre selection
    if(varPreSelectInput) {
      Xinput <- X[, attr(estKDSN, "preSelectVars"), drop=FALSE]
    }
    else{
      Xinput <- X
    }
    
    # Fit all models and save them to disk. Discard model if not used
    filePaths <- vector("character", ensembleSize)
    for(i in 1:ensembleSize) {
      
      # Bagging
      if(bagging[i]) {
        # Draw random bootstrap samples
        set.seed(seedBag[i])
        bootSamples <- sample(1:length(y), size=levelsInput*length(y), replace=TRUE)
        # Convert to list
        bootSamples <- split(bootSamples, rep(1:levelsInput, each=length(y)))
      }
      else{
        bootSamples <- NULL
      }
      
      # Draw random Fourier transformation weights
      set.seed(seedW1[i])
      posNumbers <- sample.int(.Machine$integer.max, size=levelsInput)
      set.seed(seedW2[i])
      signNumbers <- sample(c(-1, 1), size=levelsInput, replace=TRUE)
      seedrandFourInput <- signNumbers * posNumbers
      
      # Draw random dropouts
      set.seed(seedDrop1[i])
      posNumbDrop <- sample.int(.Machine$integer.max, size=levelsInput)
      set.seed(seedDrop2[i])
      signNumbersDrop <- sample(c(-1, 1), size=levelsInput, replace=TRUE)
      seedDropInput <- signNumbersDrop * posNumbDrop
      
      # Random subsets
      if(randSubsets[i]) {
        
        # Draw size of random subset
        set.seed(seedRandSubSize[i])
        randSize <- sample(1:ncol(Xinput), size=levelsInput, replace=TRUE)
        
        # Draw random covariables given size as list
        randSubInd <- lapply(1:levelsInput, function(x) {
          seedInput <- ifelse(is.null(seedRandSubDraw[i]), NULL, seedRandSubDraw[i] + x)
          set.seed(seedInput)
          sample(1:ncol(Xinput), size=randSize[x], replace=FALSE)
        }
        )
      }
      else{
        randSubInd <- NULL
      }
      
      # Fit KDSN
      ensembleModel <- fitKDSN(y = y, X = Xinput, levels = levelsInput, 
                                     Dim = DimInput, sigma = sigmaInput, 
                                     lambdaRel = lambdaRelInput, alpha = alphaInput, 
                                     info = FALSE, seedW = seedrandFourInput, 
                                     standX = standXinput, standY=standYinput,
                                     varSelect=varSelectInput, varRanking=varRankingInput,
                                     varExFreq=varExFreqInput, dropHidden=dropHiddenInput, 
                                     dropHiddenProb=dropHiddenProbInput,
                                     seedDrop=seedDropInput,
                                     baggingInd=bootSamples, randSubsetInd=randSubInd)
      
      # Include variable pre selection information, if applicable
      if(varPreSelectInput) {
        attr(ensembleModel, "preSelectVars") <- attr(estKDSN, "preSelectVars")
      }
      
      # Save bagging parameters
      attr(ensembleModel, "bagging") <- bagging[i]
      attr(ensembleModel, "seedBag") <- seedBag[i]
      # Save random subsets parameters
      attr(ensembleModel, "randSubsets") <- randSubsets[i]
      attr(ensembleModel, "seedRandSubSize") <- seedRandSubSize[i]
      attr(ensembleModel, "seedRandSubDraw") <- seedRandSubDraw[i]
      # Save random fourier transformation seeds
      attr(ensembleModel, "seedW1") <- seedW1[i]
      attr(ensembleModel, "seedW2") <- seedW2[i]
      # Save random dropout seeds
      attr(ensembleModel, "seedDrop1") <- seedDrop1[i]
      attr(ensembleModel, "seedDrop2") <- seedDrop2[i]
      
      # Save model on disk and remove it from workspace
      filePaths[i] <- paste(dirNameDisk, i, sep="")
      save(ensembleModel, file=filePaths[i], compress="xz")
      rm(ensembleModel)
      
      if(info) {
        cat("Ensemble", round(i/ensembleSize, 4)*100, "%", "\n")
      }
    }
    
    # Prepare output
    class(filePaths) <- "KDSNensembleDisk"
    attr(filePaths, "ensembleSize") <- ensembleSize
    return(filePaths)
  }
}
