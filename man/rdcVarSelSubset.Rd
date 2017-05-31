\name{rdcVarSelSubset}
\alias{rdcVarSelSubset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variable selection based on RDC with genetic algorithm (experimental)
}
\description{
Selects important variables, which have high RDC scores. A genetic algorithm is used to search the discrete space. Note that this function is still experimental.
}
\usage{
rdcVarSelSubset(x, y, k=20, s=1/6, f=sin, seedX=1:10, seedY=-c(1:10), 
rdcRep=10, popSize=100, maxiter=100, nCores=1, addInfo=TRUE)
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
  \item{rdcRep}{
    Gives the number of rdc repetitions. All repetitions are averaged per variable, to give more robust estimates. Default is to use one repetition.
  }
  \item{popSize}{
    Size of population of the genetic algorithm.
  }
  \item{maxiter}{
    Maximum number of generations to generate.
  }
  \item{nCores}{
    Number of threads used in parallelisation in \code{\link[GA]{ga}}. Default is no parallelisation.
  }
  \item{addInfo}{
    Should details of the optimization be printed? (logical scalar) Default TRUE enables        default monitoring, see \code{\link[GA]{ga}} for further details. 
  }
}
%\details{
%}
\value{
Indices of selected variables
}
\references{
David Lopez-Paz and Philipp Hennig and Bernhard Schoelkopf, (2013), 
\emph{The Randomized dependence coefficient},
Proceedings of Neural Information Processing Systems (NIPS) 26, Stateline Nevada USA, C.J.C. Burges and L. Bottou and M. Welling and Z. Ghahramani and K.Q. Weinberger (eds.) 

M. Wahde, (2008), Biological inspired methods: An introduction, WIT Press
}
\author{
Thomas Welchowski
}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{rdcPart}}, \code{\link{cancorRed}}, \code{\link{rdcSubset}}, \code{\link{rdcVarOrder}}
}
\examples{
# Generate 10 covariates
library(mvtnorm)
set.seed(3489)
X <- rmvnorm(n=200, mean=rep(0, 10))

# Generate responses based on some covariates
set.seed(-239247)
y <- 0.5*X[, 1]^3 - 2*X[, 2]^2 + X[, 3] - 1 + rnorm(200)

# Running variable selection
foundVar <- rdcVarSelSubset(x=X, y=y, seedX=1, seedY=-(1), rdcRep=1, 
popSize=80, maxiter=5)
foundVar

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
