\name{lossGCV}
\alias{lossGCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generalized cross-validation loss
}
\description{
Calculates the generalized cross-validation score given number of observations, model deviance, effective degrees of freedom and gamma. 
}
\usage{
lossGCV(n, Dev, trA, gammaPar=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of observations (integer scalar).
}
  \item{Dev}{
Deviance of the kernel ridge regression (numeric scalar).
}
  \item{trA}{
Effective degrees of freedom (numeric scalar).
}
  \item{gammaPar}{
Weighting parameter (numeric scalar), which specifies how the generalized cross-validation score is penalized by the effective degrees of freedom. Default value is 1.
}
}
%\details{
%
%}
\value{
Generalized cross-validation loss (numeric scalar).
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
\code{\link{lossApprox}}, \code{\link{lossCvKDSN}}, \code{\link{lossSharedCvKDSN}}, 
\code{\link{lossSharedTestKDSN}}
}

\examples{
# Simple example based on help pages of mgcv::gam
# GCV Score is the same as used in function mgcv::mgcv
library(mgcv)
dat <- gamSim(1, n=100, dist="normal", scale=2)
gamMod <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),data=dat)
calcGCV <- lossGCV(n=100, Dev=deviance(gamMod), trA=sum(gamMod$edf), gammaPar=1)
all.equal(calcGCV, gamMod$gcv.ubre) # TRUE
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line