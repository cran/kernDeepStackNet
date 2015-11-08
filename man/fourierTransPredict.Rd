\name{fourierTransPredict}
\alias{fourierTransPredict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction based on random Fourier transformation
}
\description{
Based on prior given weights, the estimated Fourier transformation is applied to new data.
}
\usage{
fourierTransPredict(newx, rW)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{newx}{
New data design matrix.
}
  \item{rW}{
Prior drawn random weight matrix.
}
}
%\details{
%}
\value{
Numeric transformed data matrix with dimension \code{2*Dim x n}.
}
\references{
Po-Seng Huang and Li Deng and Mark Hasegawa-Johnson and Xiaodong He, (2013), 
\emph{Random Features for kernel deep convex network},
Proceedings IEEE International Conference on Acoustics, Speech, and 
Signal Processing (ICASSP)
}
\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}
}

\note{
This function is not intended to be called directly by the user. Should only be used by experienced users, who want to customize the model. It is called in the estimation process of the kernel deep stacking network, e.g. \code{\link{fitKDSN}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{randomFourierTrans}}
}

\examples{
# Generate data matrix
set.seed(50)
X <- matrix(rnorm(100*3), ncol=3)

# Apply a random Fourier transformation of higher dimension
rft <- randomFourierTrans(X=X, Dim=3, sigma=1, seedW=0)

# New data matrix
set.seed(100)
Xnew <- matrix(rnorm(100*3), ncol=3)

# Apply same Fourier transformation on new data
newZ <- fourierTransPredict(newx=Xnew, rW=rft$rW)
head(newZ)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line