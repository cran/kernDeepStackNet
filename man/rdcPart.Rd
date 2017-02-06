\name{rdcPart}
\alias{rdcPart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Randomized dependence coefficient partial calculation
}
\description{
  Calculates the randomized dependence coefficient based on interim results. It is a generalized dependence measure based on maximum correlation of random non-linear projections.
}
\usage{
  rdcPart(subsetX, xTrans, yTrans, s=1/6, f=sin, randX)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{subsetX}{
    Subset of the covariate matrix as indices (integer vector).
  }
  \item{xTrans}{
    Transformed matrix to the [0, 1] scale (numeric matrix).
  }
  \item{yTrans}{
   Random, non-linear projection of the response (numeric vector).
  }
  \item{s}{
    Variance of the random weights. Default is 1/6.
  }
   \item{f}{
    Non-linear transformation function. Default is \code{\link[base]{sin}}.
  }
   \item{randX}{
    Random weights (numeric vector).
  }
}
\details{
  This function allows for more efficient calculation than the complete calculation by excluding repetitive calculations.
}

\references{
David Lopez-Paz et. al, (2013), The randomized dependence coefficient, Proceedings of Advances in Neural Information Processing Systems 26 (NIPS)
}

\value{
  Value of randomized dependence coefficient (numeric scalar).
}

\author{
  Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}
}

\note{
  This function is a help function within variable selection. It is given for experts for model customization. It is recommended to use instead function \code{\link{rdcVarOrder}.
}
}

\seealso{
\code{\link{rdcVarOrder}}, \code{\link{cancorRed}}
}

%\examples{
 %%
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
