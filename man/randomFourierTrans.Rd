\name{randomFourierTrans}
\alias{randomFourierTrans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Random Fourier transformation
}
\description{
Calculates the random Fourier transformation using a Gaussian kernel, given the original data.
}
\usage{
randomFourierTrans(X, Dim, sigma, seedW=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Original data design matrix. All factors have to be encoded, e.g. dummy coding.
}
  \item{Dim}{
Specifies the dimension of the random Fourier transformation (integer scalar).
}
  \item{sigma}{
Variance of the Gaussian kernel (positive numeric scalar).
}
  \item{seedW}{
Random seed for drawing from the multivariate normal distribution (integer scalar).
}
}
\details{
First a random weight matrix is drawn from the multivariate normal distribution. Then the Data is linear transformed. The linear transformed data is mapped nonlinear by applying cosinus and sinus functions. The matrix multiplaction t(Z) Z approximates the Gaussian radial basis function kernel matrix. The dimension of t(Z) Z is always n x n. The higher the dimension argument Dim, the more accurate the results.
}
\value{
Numeric transformed data matrix with dimension 2*Dim x n.
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
This function is not intended to be called directly by the user. Should only be used by experienced users, who want to customize the model. It is called in the estimation process of the kernel deep stacking network, e. g. \code{\link{fitKDSN}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fourierTransPredict}}
}

\examples{
# Generate data matrix
X <- data.frame(rnorm(100), rnorm(100), rnorm(100), rnorm(100), rnorm(100), 
factor(sample(c("a", "b", "c", "d", "e"), 100, replace=TRUE)))
X <- model.matrix(object=~., data=X) 
# Exclude intercept
X <- X[, -1]

# Apply a random Fourier transformation of lower dimension
rft <- randomFourierTrans(X=X, Dim=2, sigma=1, seedW=0)

# Transformed data
rft$Z

# Used weight matrix
rft$rW
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line