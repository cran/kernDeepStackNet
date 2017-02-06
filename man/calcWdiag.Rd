\name{calcWdiag}
\alias{calcWdiag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of weight matrix
}
\description{
Calculates the weight matrix in the pseudo iterative reweighted least squares step of a generalized linear model. It is based on the variance function and the first derivative of the link function. The values are calculated at the actual expected values. The dimension of the matrix is n x n (n = number of observations). For details see reference.
}
\usage{
calcWdiag(varMu, gDerivMu)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{varMu}{
    Value of the variance function of the exponential family on the expected value mu (numeric scalar).
  }
  \item{gDerivMu}{
    Value of the first derivative of the link function evaluated at the expected value. (numeric scalar).
  }
}
%\details{
%
%}
\value{
Weight matrix (numeric vector).
}
\references{
Simon N. Wood, (2006), 
\emph{Generalized Additive Models: An Introduction with R},
Taylor \& Francis Group LLC
}
\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}
}

%\note{
%%This is a help function to calculate the loss function of the kernel deep stacking network, e.g. \code{\link{lossKDSN}}. It is not intended to be called directly by the user, but is given for experienced users to customize the model.
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{calcTrA}}, \code{\link{varMu}}, \code{\link{gDerivMu}}
}

%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line