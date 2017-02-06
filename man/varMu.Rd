\name{varMu}
\alias{varMu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variance function evaluated at expected value
}
\description{
Evaluates the variance function at the expected value. Only implemented for Gaussian distribution.
}
\usage{
varMu(mu)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{mu}{
Expected values conditional on the design matrix. These are the fitted values of the model (numeric vector).
}
}
%\details{
%}
\value{
Variance function of exponential family evaluated at the fitted values (numeric vector).
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
\code{\link{calcTrA}}, \code{\link{calcWdiag}}, \code{\link{gDerivMu}}
}

%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line