\name{rdcVarOrder}
\alias{rdcVarOrder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Variable ordering using randomized dependence coefficients (experimental)
}
\description{
  Variable selection for KDSN. All variables are univariately compared to the responses with the randomized dependence coefficient (RDC). The indices of the variables are ordered accordingly. Variables with low RDC may be discarded.
}
\usage{
  rdcVarOrder(x, y, k=20, s=1/6, f=sin, seedX=NULL, 
            seedY=NULL, nCores=1, info=FALSE, cutoff=0, rdcRep=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Covariates data (numeric matrix).
  }
  \item{y}{
    Responses (numeric matrix).
  }
  \item{k}{
    Number of random features (integer scalar).
  }
  \item{s}{
    Variance of the random weights. Default is 1/6.
  }
  \item{f}{
    Non-linear transformation function. Default is \code{\link[base]{sin}}.
  }
  \item{seedX}{
    Random number seed of normal distributed weights for covariates (integer scalar). Default is to randomly draw weights.
  }
  \item{seedY}{
     Random number seed of normal distributed weights for responses (integer scalar). Default is to randomly draw weights.
  }
  \item{nCores}{
    Number of threads used. If greater than one, parallel processing using the parallel package is used.
  }
  \item{info}{
    Should additional infos on the progress of the function be given? (logical scalar) Default is not to give additional information.
  }
  \item{cutoff}{
    Gives the empirical cut-off value (numeric scalar). Variables below this threshold will be discarded. Default is to include all variables.
  }
  \item{rdcRep}{
    Gives the number of rdc repetitions. All repetitions are averaged per variable, to give more robust estimates. Default is to use one repetition.
  }
}

\details{
Covariates are ranked according to their dependence with the response variable. Note that this function is still experimental.
}

\references{
David Lopez-Paz and Philipp Hennig and Bernhard Schoelkopf, (2013), 
\emph{The Randomized dependence coefficient},
Proceedings of Neural Information Processing Systems (NIPS) 26, Stateline Nevada USA, C.J.C. Burges and L. Bottou and M. Welling and Z. Ghahramani and K.Q. Weinberger (eds.) 
}

\value{
  Ordered indices of original covariates in increasing order of importance. The first variable is the least important. 
}

\author{
  Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}
}

%\note{
%
%}
  
\seealso{
  \code{\link{rdcPart}}, \code{\link{cancorRed}}
}
  
\examples{
#############################
# Cubic noisy association

# Generate 10 covariates
library(mvtnorm)
set.seed(3489)
X <- rmvnorm(n=200, mean=rep(0, 10))

# Generate responses based on some covariates
set.seed(-239247)
y <- 0.5*X[, 1]^3 - 2*X[, 2]^2 + X[, 3] - 1 + rnorm(200)

# Variable selection with RDC
selectedInd <- rdcVarOrder(x=X, y=y, seedX=1, seedY=2, cutoff=0.7)
selectedInd
# -> If the numbers of variables should be reduced from 10 to 3, 
# then all important variables are selected.

# With more repetitions and different random transformations
selectedInd <- rdcVarOrder(x=X, y=y, seedX=1:25, seedY=-(1:25), cutoff=0.7, rdcRep=25)
selectedInd
# -> Gives identical result as one repetition
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
  