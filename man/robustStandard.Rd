\name{robustStandard}
\alias{robustStandard}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robust standardization
}
\description{
Scales the data matrix with the median and median absolute deviation.
}
\usage{
robustStandard(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{X}{
Data design matrix.
}
}
%\details{
%}
\value{
Numeric transformed data matrix with same dimension as original data.
}

\references{
Ricardo A. Maronna and R. Douglas Martin and Victor J. Yohai, (2006),
\emph{Robust Statistics: Theory and Methods},
John Wiley \& Sons, Ltd.
}

%@book{robustStat,
% title={Robust Statistics: Theory and Methods},
% author={Ricardo A. Maronna and R. Douglas Martin and Victor J. Yohai},
% year={2006},
% publisher={John Wiley \& Sons, Ltd.}
%}

\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}
}

%\note{
%This function is not intended to be called directly by the user. Should only be used by %experienced users, who want to customize the model. It is called in the estimation process %of the kernel deep stacking network, e.g. \code{\link{fit_k_DSN_rft}}.
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%\code{\link{fourierTransPredict}}
%}

\examples{
# Generate data matrix
set.seed(150)
X <- matrix(rnorm(100*3), ncol=3)

# Robust standardization
scaledX <- robustStandard(X=X)

# Median equals 0
all.equal(sapply(1:3, function(x) median(scaledX[,x])), rep(0, 3))

# MAD equals 1
all.equal(sapply(1:3, function(x) mad(scaledX[,x], constant=1)), rep(1, 3))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line