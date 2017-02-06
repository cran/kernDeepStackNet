\name{calcTrA}
\alias{calcTrA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the trace of the hat matrix
}
\description{
The trace of the hat matrix corresponds to the effective degrees of freedom (edf) of a generalized additive model. The edf is used as an internal measure of model complexity.
}
\usage{
calcTrA(W, X, lambda=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{W}{
Weight matrix of the pseudo iterated least squares algorithm. See function \code{\link{calcWdiag}}
}
\item{X}{
Design matrix of the covariates.
}
\item{lambda}{
Regularization parameter of kernel ridge regression.
}
}
%\details{
%
%}
\value{
Effective degrees of freedom of a generalized additive model with regularization (numeric scalar).
}
\references{
Simon N. Wood, (2006), 
\emph{Generalized Additive Models: An Introduction with R},
Taylor \& Francis Group LLC
}
\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}
}

\note{
This function is not intended to be called directly by the user. Should only be used by experienced users, who want to customize the model. %%It is called in the model selection process of the kernel deep stacking network, e.g. \code{\link{tuneKDSN}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{varMu}}, \code{\link{calcWdiag}}, \code{\link{gDerivMu}}
}

%\examples{
%
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line