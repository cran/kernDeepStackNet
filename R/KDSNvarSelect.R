# Canonical correlation analysis
cancorRed <- function(x, y, xcenter = TRUE, ycenter = TRUE)  {
#   x <- as.matrix(x)
#   y <- as.matrix(y)
  if ((nr <- nrow(x)) != nrow(y)) 
    stop("unequal number of rows in 'cancor'")
  ncx <- ncol(x)
  ncy <- ncol(y)
  if (!nr || !ncx || !ncy) 
    stop("dimension 0 in 'x' or 'y'")
  if (is.logical(xcenter)) {
    if (xcenter) {
      xcenter <- colMeans(x)
      x <- x - rep(xcenter, rep.int(nr, ncx))
    }
    else xcenter <- rep.int(0, ncx)
  }
  else {
    xcenter <- rep_len(xcenter, ncx)
    x <- x - rep(xcenter, rep.int(nr, ncx))
  }
  if (is.logical(ycenter)) {
    if (ycenter) {
      ycenter <- colMeans(y)
      y <- y - rep(ycenter, rep.int(nr, ncy))
    }
    else ycenter <- rep.int(0, ncy)
  }
  else {
    ycenter <- rep_len(ycenter, ncy)
    y <- y - rep(ycenter, rep.int(nr, ncy))
  }
  qx <- qr(x)
  qy <- qr(y)
  dx <- qx$rank
  if (!dx) 
    stop("'x' has rank 0")
  dy <- qy$rank
  if (!dy) 
    stop("'y' has rank 0")
  z <- svd(qr.qty(qx, qr.qy(qy, diag(1, nr, dy)))[1L:dx, , 
                                                  drop = FALSE], dx, dy)
  return(z$d[1])
}

# RDC evaluation of covariates
rdcPart <- function(subsetX, xTrans, yTrans, s=1/6, f=sin, randX) {
  xTrans <- xTrans[, c(subsetX, ncol(xTrans))]
  xTrans <- s/ncol(xTrans) * 
    xTrans %*% matrix(randX, ncol(xTrans))
  return(cancorRed(cbind(f(xTrans), 1), cbind(yTrans, 1)))
}

# Function for variable selection based on rdc
rdcVarOrder <- function(x, y, k=20, s=1/6, f=sin, seedX=NULL, 
                      seedY=NULL, nCores=1, 
                      info=FALSE, cutoff=0, rdcRep=1) {
  
  if(rdcRep==1) {
  
    # Transformation to scale [0, 1]
    # x <- cbind(apply(as.matrix(x),2,function(u)rank(u)/length(u)),1)
    # y <- cbind(apply(as.matrix(y),2,function(u)rank(u)/length(u)),1)
    x <- colRanks(as.matrix(x), preserveShape = TRUE, ties.method="average")
    x <- cbind(x / nrow(x), 1)
    y <- colRanks(as.matrix(y), preserveShape = TRUE, ties.method="average")
    y <- cbind(y / nrow(y), 1)
    
    # Draw random numbers
    set.seed(seedX)
    xRand <- rnorm(2*k)
    set.seed(seedY)
    yRand <- rnorm(ncol(y)*k)
    
    # Calculate random nonlinear projection of responses
    y <- s/ncol(y) * y %*% matrix(yRand, ncol(y))
    y <- f(y)
  
    # Create indices design
    rdcDesign <- 1:(dim(x)[2]-1)
    
    # Calculation of rdc scores
    if(nCores==1) {
      rdcScores <- vector("numeric", length(rdcDesign))
      for(i in 1:length(rdcDesign)) {
        rdcScores[i] <- rdcPart (subsetX=rdcDesign[i], xTrans=x, yTrans=y, s=s, f=f, 
                                 randX=xRand)
        if(info) {cat("Variables:", round(i/length(rdcDesign), 4)*100, "% done", "\n")}
      }
    }
    else{
      socketClust <- makeCluster(nCores) 
      clusterExport(cl = socketClust, varlist=c("rdcPart", "cancorRed",
                                                "x", "y", "s", "f", "xRand", 
                                                "rdcDesign"), envir=environment())
      rdcScores <- parSapply(cl=socketClust, X=1:length(rdcDesign), 
                             FUN=function(i) 
                               rdcPart (subsetX=rdcDesign[i], xTrans=x, yTrans=y, 
                                        s=s, f=f, randX=xRand))
      stopCluster(socketClust)
    }
    
    # Select variables with higher values than the empirical cut-off quantile
    # Order is increasing -> First exclude the beginning index and continue
    checkRes <- rdcScores >= quantile(rdcScores, prob = cutoff)
    output <- which(checkRes) [order( rdcScores[checkRes] ) ]
    return(output)
  }
  
  # With repetitions
  else{

    # Transformation to scale [0, 1]
    # x <- cbind(apply(as.matrix(x),2,function(u)rank(u)/length(u)), 1)
    # yPre <- cbind(apply(as.matrix(y),2,function(u)rank(u)/length(u)), 1)
    x <- colRanks(as.matrix(x), preserveShape = TRUE, ties.method="average")
    x <- cbind(x / nrow(x), 1)
    yPre <- colRanks(as.matrix(y), preserveShape = TRUE, ties.method="average")
    yPre <- cbind(yPre / nrow(yPre), 1)
    
    # Create indices design
    rdcDesign <- 1:(dim(x)[2]-1)
    
    # Calculation of rdc scores
    if(nCores==1) {
      rdcScores <- matrix(NA, nrow=length(rdcDesign), ncol=rdcRep)
      
      for(j in 1:rdcRep) {
        # Draw random numbers
        set.seed(seedX[j])
        xRand <- rnorm(2*k)
        set.seed(seedY[j])
        yRand <- rnorm(ncol(yPre)*k)
        
        # Calculate random nonlinear projection of responses
        y <- s/ncol(yPre) * yPre %*% matrix(yRand, ncol(yPre))
        y <- f(y)
      
        for(i in 1:length(rdcDesign)) {
            rdcScores[i, j] <- rdcPart (subsetX=rdcDesign[i], xTrans=x, yTrans=y, s=s, f=f, randX=xRand)
            if(info) {cat("Variables:", round(i/length(rdcDesign), 4)*100, "% done", "\n")}
          }
        if(info) {cat("Repetitions:", round(j/rdcRep, 4)*100, "% done", "\n")}
      }
      
      # Average over all repetitions
      rdcScores <- rowMeans(rdcScores)
    }
    else{
      # Define function to parallelize
      tempFun <- function(j) {
        rdcScores <- vector("numeric", length(rdcDesign))
          # Draw random numbers
          set.seed(seedX[j])
          xRand <- rnorm(2*k)
          set.seed(seedY[j])
          yRand <- rnorm(ncol(yPre)*k)
          
          # Calculate random nonlinear projection of responses
          y <- s/ncol(yPre) * yPre %*% matrix(yRand, ncol(yPre))
          y <- f(y)
          
          for(i in 1:length(rdcDesign)) {
            rdcScores[i] <- rdcPart (subsetX=rdcDesign[i], xTrans=x, 
                                        yTrans=y, s=s, f=f, randX=xRand)
          }
        return(rdcScores)
      }
      
      socketClust <- makeCluster(nCores) 
      clusterExport(cl = socketClust, varlist=c("tempFun", "rdcPart", "cancorRed",
                                                "y", "x", "s", "f", "rdcDesign",
                                                "seedX", "seedY"), 
                    envir=environment())
      rdcScoresPre <- parSapply(cl=socketClust, X=1:rdcRep, FUN=tempFun)
      stopCluster(socketClust)
      rdcScores <- rowMeans(rdcScoresPre)
    }
    
    # Select variables with higher values than the empirical cut-off quantile
    # Order is increasing -> First exclude the beginning index and continue
    checkRes <- rdcScores >= quantile(rdcScores, prob = cutoff)
    output <- which(checkRes) [order( rdcScores[checkRes] ) ]
    return(output)
  }
  
}

# Function for variable pre selection
rdcSubset <- function(binCode, x, y, k=20, s=1/6, f=sin, seedX=NULL, 
                      seedY=NULL, rdcRep=1, trans0to1=TRUE) {
  # BinCode must include at least one covariate, otherwise the measure is defined as zero
  if(sum(binCode)==0) {
    return(0)
  }
  
  # Conversion of input
  subsetX <- which(binCode==1)
  
  if(trans0to1) {
    # Transformation to scale [0, 1]
    # x <- cbind(apply(as.matrix(x),2,function(u)rank(u)/length(u)), 1)
    # y <- cbind(apply(as.matrix(y),2,function(u)rank(u)/length(u)), 1)
    x <- colRanks(as.matrix(x), preserveShape = TRUE, ties.method="average")
    x <- cbind(x / nrow(x), 1)
    yPre <- colRanks(as.matrix(y), preserveShape = TRUE, ties.method="average")
    yPre <- cbind(yPre / nrow(yPre), 1)
  }
  else{
    yPre <- y
  }

  # Calculation of rdc scores
  rdcScores <- vector("numeric", length=rdcRep)
  for(j in 1:rdcRep) {
    
    # Draw random numbers
    set.seed(seedX[j])
    xRand <- rnorm((length(subsetX)+1)*k)
    set.seed(seedY[j])
    yRand <- rnorm(ncol(yPre)*k)
    
    # Calculate random nonlinear projection of responses
    y <- s/ncol(yPre) * yPre %*% matrix(yRand, ncol(yPre))
    y <- f(y)
    
    # Calculate rdc scores
    rdcScores[j] <- rdcPart (subsetX=subsetX, xTrans=x, yTrans=y, s=s, f=f, randX=xRand)
  }
  
  # Average over all repetitions
  rdcScore <- sum(rdcScores) / length(rdcScores)
  return(rdcScore)
}

rdcVarSelSubset <- function(x, y, k=20, s=1/6, f=sin, seedX=1:10, 
                            seedY=-c(1:10), rdcRep=10, 
                            popSize=100, maxiter=100, nCores=1, addInfo=TRUE) {
  
  # Transform to [0, 1]
  x <- colRanks(as.matrix(x), preserveShape = TRUE, ties.method="average")
  x <- cbind(x / nrow(x), 1)
  y <- colRanks(as.matrix(y), preserveShape = TRUE, ties.method="average")
  y <- cbind(y / nrow(y), 1)
  
#   # Definition of function to optimize
#   gainFunc <- function(i, x, y, rdcRep, seedX, seedY) {
# 
# #    if(sum(i)!=0) {
# #    }
# #    else{
# #      0
# #    }
# 
#     rdcSubset(binCode=i, x=x, y=y, rdcRep=rdcRep, seedX=seedX, seedY=seedY, 
#               trans0to1=FALSE)
#   }
  if(nCores==1) {
  
    if(addInfo){
      gaRes <- ga(type="binary", fitness=rdcSubset, x=x, y=y, rdcRep=rdcRep, 
                  seedX=seedX, seedY=seedY, trans0to1=FALSE,
                  nBits=ncol(x)-1, parallel=FALSE, 
                  popSize=popSize, maxiter=maxiter)
    }
    else{
      gaRes <- ga(type="binary", fitness=rdcSubset, x=x, y=y, rdcRep=rdcRep, 
                  seedX=seedX, seedY=seedY, trans0to1=FALSE,
                  nBits=ncol(x)-1, parallel=FALSE, 
                  popSize=popSize, maxiter=maxiter, monitor=FALSE)
    }
  }
  else{
    socketClust <- makeCluster(nCores) 
    clusterExport(cl = socketClust, varlist=c("rdcSubset", "rdcPart", "cancorRed"), envir=environment())
    gaRes <- ga(type="binary", fitness=rdcSubset, x=x, y=y, rdcRep=rdcRep, 
                seedX=seedX, seedY=seedY, trans0to1=FALSE,
                nBits=ncol(x)-1, parallel=TRUE, 
                popSize=popSize, maxiter=maxiter)
    stopCluster(socketClust)
  }
  return(which(gaRes@solution[1, ]==1))
}
