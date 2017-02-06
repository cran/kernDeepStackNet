\name{rdcSubset}
\alias{rdcSubset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Randomized dependence coefficients score on given subset
}
\description{
  Variable pre selection scoring for KDSN. Estimates the RDC score for a subset of variables. 
}
\usage{
  rdcSubset(binCode, x, y, k=20, s=1/6, f=sin, seedX=NULL, seedY=NULL, 
  rdcRep=1, trans0to1=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{binCode}{
    Specifies which set of variables of the covariates is used to explain the responses (binary vector). One to assiged inclusion and zero excludes variables.
  }
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
  \item{rdcRep}{
    Gives the number of rdc repetitions. All repetitions are averaged per variable, to give more robust estimates. Default is to use one repetition.
  }
    \item{trans0to1}{
    Should the design matrix and response be transformed to the interval [0, 1]? (Logical). If the data is available in this for form, it can be evaluated much faster.
  }
}

\details{
Covariates are ranked according to their dependence with the response variable. 
}

\references{
David Lopez-Paz and Philipp Hennig and Bernhard Schoelkopf, (2013), 
\emph{The Randomized dependence coefficient},
Proceedings of Neural Information Processing Systems (NIPS) 26, Stateline Nevada USA, C.J.C. Burges and L. Bottou and M. Welling and Z. Ghahramani and K.Q. Weinberger (eds.) 
}

\value{
RDC score (numeric scalar).
}

\author{
  Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}
}

%\note{
%
%}
  
\seealso{
  \code{\link{rdcPart}}, \code{\link{cancorRed}}, \code{\link{rdcVarOrder}}, \code{\link{rdcVarSelSubset}}
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

# Score of true subset
scoreTrue <- rdcSubset(binCode=c(rep(1, 3), rep(0, 7)), 
x=X, y=y, seedX=1:10, seedY=-(1:10), rdcRep=10)
scoreTrue

# Only unneccessary variables
scoreFalse <- rdcSubset(binCode=c(rep(0, 3), rep(1, 7)), 
x=X, y=y, seedX=1:10, seedY=-(1:10), rdcRep=10)
scoreFalse

# One important two important variables and some non causal variables
scoreMix <- rdcSubset(binCode=c(1, 0, 1, rep(0, 3), rep(1, 4)), 
x=X, y=y, seedX=1:10, seedY=-(1:10), rdcRep=10)
scoreMix

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
  