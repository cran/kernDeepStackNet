# General sequential one dimensional optimization of f(x), x \in R^d, f(x) \in R
optimize1dMulti <- function (f_input, lower, upper, maxRuns=3, repetitions=5, 
                             tol_input=.Machine$double.eps^0.25, x_0=NULL, addInfo=TRUE,
                             nCores=1, envir=parent.frame(), directUse=TRUE, OptTypePar="") {
  
  # Check if first argument of function is x
  stopifnot(formalArgs (f_input) [1]=="x")

  # Rerun optimization with different starting values
  Results <- vector("list", repetitions)
  # dimension <- dim(interval_matrix) [2]
  dimension <- length(lower)
  if(nCores==1) {
    for(j in 1:repetitions) {
      
      if(j > 1 | is.null(x_0)) {
        # Set initial starting value: Random vector x nid ~ U (a_i, b_i)
        x_0 <- sapply(1:dimension, function (x) runif(n=1, min=lower[x], max=upper[x]))
      }
      x_0_alt <- x_0
      liste_vektoren <- vector("list", dimension)
      
      if(dimension==1) {
        liste_vektoren [[1]] <- expression(x)
      }
      
      if(dimension==2) {
        liste_vektoren [[1]] <- expression(c(x, x_0 [2]))
        liste_vektoren [[2]] <- expression(c(x_0 [1], x))
      }
      
      if(dimension>=3) {
        liste_vektoren [[1]] <- expression(c(x, x_0 [2:dimension]))
        liste_vektoren [[dimension]] <- expression(c(x_0 [1:(dimension-1)], x))
        for(i in 1:(dimension-2)) {
          liste_vektoren [[i+1]] <- substitute(c(x_0 [1:i], x, x_0 [(2+i):dimension]),list(i=i))
        }
      }
      
      # Univariate optimization over one variable given all other variables
      i <- 1
      whileCondition <- TRUE
      stepRuns <- 0
      f_input_x_0 <- f_input(x_0)
      f_input_x_0_alt <- f_input_x_0
      x_0_alt <- x_0
      
      while (whileCondition) {
        # Conditional optimization
        for(i in 1:dimension) {
          # erg <- optimize(f=function (x) f_input( eval( liste_vektoren [[i]] ) ), interval=interval_matrix [, i], tol=tol_input)
          erg <- optimize(f=function (x) f_input( eval( liste_vektoren [[i]] ) ), interval=c(lower[i], upper[i]), tol=tol_input)
          x_0 [i] <- erg$minimum
          if(addInfo) {
            cat("optimize1dMulti", "parameter", i, "\n")
          }
        }
        
        # Condition for stopping
        stepRuns <- stepRuns + 1
        f_input_x_0 <- erg$objective
        if(stepRuns == maxRuns) {
          whileCondition <- FALSE
        }
        else {
          if(abs(f_input_x_0_alt) != 0 & sum(abs(x_0_alt)) != 0) {
            whileCondition <- (abs(f_input_x_0 - f_input_x_0_alt) / abs(f_input_x_0_alt) >= tol_input) & (sum(abs(x_0 - x_0_alt)) / sum(abs(x_0_alt)) >= tol_input)
          }
          else {
            if(abs(f_input_x_0_alt) == 0 & sum(abs(x_0_alt)) != 0) {
              whileCondition <- (abs(f_input_x_0 - f_input_x_0_alt) >= tol_input) & (sum(abs(x_0 - x_0_alt)) / sum(abs(x_0_alt)) >= tol_input)
            }
            if(abs(f_input_x_0_alt) != 0 & sum(abs(x_0_alt)) == 0) {
              whileCondition <- (abs(f_input_x_0 - f_input_x_0_alt) / abs(f_input_x_0_alt) >= tol_input) & (sum(abs(x_0 - x_0_alt)) >= tol_input)
            }
            if(abs(f_input_x_0_alt) == 0 & sum(abs(x_0_alt)) == 0) {
              whileCondition <- (abs(f_input_x_0 - f_input_x_0_alt) >= tol_input) & (sum(abs(x_0 - x_0_alt)) >= tol_input)
            }
          }
        }
        f_input_x_0_alt <- f_input_x_0
        x_0_alt <- x_0
        if(addInfo) {
          cat("optimize1dMulti", "run", stepRuns, "\n")
        }
      }
      Results [[j]] <- list(minimum=x_0, objective=f_input_x_0)
      if(addInfo) {
        cat("optimize1dMulti", "repetition", j, "\n")
      }
    }
  }
  else{

    # Help function for parallelisation
    tempFunc <- function(j) {
      
      if(j > 1 | is.null(x_0)) {
        # Set initial starting value: Random vector x nid ~ U (a_i, b_i)
        x_0 <- sapply(1:dimension, function (x) runif(n=1, min=lower[x], max=upper[x]))
      }
      x_0_alt <- x_0
      liste_vektoren <- vector("list", dimension)
      
      if(dimension==1) {
        liste_vektoren [[1]] <- expression(x)
      }
      
      if(dimension==2) {
        liste_vektoren [[1]] <- expression(c(x, x_0 [2]))
        liste_vektoren [[2]] <- expression(c(x_0 [1], x))
      }
      
      if(dimension>=3) {
        liste_vektoren [[1]] <- expression(c(x, x_0 [2:dimension]))
        liste_vektoren [[dimension]] <- expression(c(x_0 [1:(dimension-1)], x))
        for(i in 1:(dimension-2)) {
          liste_vektoren [[i+1]] <- substitute(c(x_0 [1:i], x, x_0 [(2+i):dimension]),list(i=i))
        }
      }

      # Univariate optimization of one variable given all other variables
      i <- 1
      whileCondition <- TRUE
      stepRuns <- 0
      f_input_x_0 <- f_input(x_0)
      f_input_x_0_alt <- f_input_x_0
      x_0_alt <- x_0
      
      while (whileCondition) {
        # Conditional optimization
        for(i in 1:dimension) {
#          erg <- optimize(f=function (x) f_input( eval( liste_vektoren [[i]] ) ), interval=interval_matrix [, i], tol=tol_input)
          erg <- optimize(f=function (x) f_input( eval( liste_vektoren [[i]] ) ), interval=c(lower[i], upper[i]), tol=tol_input)
          x_0 [i] <- erg$minimum
          if(addInfo) {
            cat("optimize1dMulti", "parameter", i, "\n")
          }
        }
        
        # Condition for stopping
        stepRuns <- stepRuns + 1
        f_input_x_0 <- erg$objective
        if(stepRuns == maxRuns) {
          whileCondition <- FALSE
        }
        else {
          if(abs(f_input_x_0_alt) != 0 & sum(abs(x_0_alt)) != 0) {
            whileCondition <- (abs(f_input_x_0 - f_input_x_0_alt) / abs(f_input_x_0_alt) >= tol_input) & (sum(abs(x_0 - x_0_alt)) / sum(abs(x_0_alt)) >= tol_input)
          }
          else {
            if(abs(f_input_x_0_alt) == 0 & sum(abs(x_0_alt)) != 0) {
              whileCondition <- (abs(f_input_x_0 - f_input_x_0_alt) >= tol_input) & (sum(abs(x_0 - x_0_alt)) / sum(abs(x_0_alt)) >= tol_input)
            }
            if(abs(f_input_x_0_alt) != 0 & sum(abs(x_0_alt)) == 0) {
              whileCondition <- (abs(f_input_x_0 - f_input_x_0_alt) / abs(f_input_x_0_alt) >= tol_input) & (sum(abs(x_0 - x_0_alt)) >= tol_input)
            }
            if(abs(f_input_x_0_alt) == 0 & sum(abs(x_0_alt)) == 0) {
              whileCondition <- (abs(f_input_x_0 - f_input_x_0_alt) >= tol_input) & (sum(abs(x_0 - x_0_alt)) >= tol_input)
            }
          }
        }
        f_input_x_0_alt <- f_input_x_0
        x_0_alt <- x_0
        if(addInfo) {
          cat("optimize1dMulti", "run", stepRuns, "\n")
        }
      }
      Results1 <- list(minimum=x_0, objective=f_input_x_0)
      if(addInfo) {
        cat("optimize1dMulti", "repetition", j, "\n")
      }
      return(Results1)
    }
    
    if(!directUse){
      
      # Subroutine to direct minimization
      if(OptTypePar=="tuneKDSN" ||  OptTypePar=="tuneLevelKDSN") {
        localEnvir <- environment()
        clusterExport(cl=envir$cl, varlist=ls(), envir=localEnvir)
        clusterExport(cl = envir$cl, varlist = objects(envir), envir=envir)
        Results <- parLapply(cl = envir$cl, X=1:repetitions, fun=tempFunc)
      }
      
      # Subroutine to Mbo minimization
      if(OptTypePar=="mbo1d") {
        localEnvir <- environment()
        clusterExport(cl=envir$envir$cl, varlist=ls(), envir=localEnvir)
        Results <- parLapply(cl = envir$envir$cl, X=1:repetitions, fun=tempFunc)
      }
    }
    else{
      cl <- makeCluster(nCores)
      localEnvir <- environment()
      clusterExport(cl=cl, varlist=ls(all.names=TRUE), envir=localEnvir)
      Results <- parLapply(cl = cl, X=1:repetitions, fun=tempFunc)
      stopCluster(cl=cl)
    }
  }
    
  # Choose best iteration among the random starting values
  Index <- which.min(sapply(1:repetitions, function (x) Results [[x]]$objective))
  Output <- list(minimum=Results [[Index]]$minimum, objective=Results [[Index]]$objective)
  return(Output)
}
