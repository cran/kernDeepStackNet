\name{lossKDSNgivenModel}
\alias{lossKDSNgivenModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kernel deep stacking network loss function based on model fit
}
\description{
Computes the generalized cross-validation (gcv) score of a tuning parameter configuration given the data. The tuning parameters consists of dimension of the random Fourier transform (=D), the variance of the Gaussian kernels (=sigma) and the regularization parameter of the kernel ridge regression (=lambda). It is a wrapper around the function \code{\link{lossGCV}}. In contrast to the function \code{\link{lossKDSN}}, there is no need to refit the model.
}
\usage{
lossKDSNgivenModel (KDSNfit, y, X, gammaPar=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{KDSNfit}{
Object of class "KDSN" \code{\link{fitKDSN}}.
}
  \item{y}{
Numeric response vector of the outcome. Should be formated as a one column matrix.
}
  \item{X}{
Numeric design matrix of the covariates. All factors have to be prior encoded.
}
  \item{gammaPar}{
Weighting parameter (numeric scalar), which specifies how the generalized cross-validation score is penalized by the effective degrees of freedom. Default value is 1.
}
}
%\details{
%
%}
\value{
Numeric scalar with gives the generalized cross-validation score. The kernel deep stacking network used to calculate the score is available as attribute.
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
This function is not intended to be called directly by the user. Should only be used by experienced users, who want to customize the model. It is called in the model selection process of the kernel deep stacking network, e.g. \code{\link{tuneKDSN}}. This version is faster than \code{\link{lossKDSN}}, because it does not require a refit.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lossKDSN}}, \code{\link{lossCvKDSN}}, \code{\link{lossGCV}}
}

\examples{
####################################
# Example with simple binary outcome

# Generate covariate matrix
sampleSize <- 100
X <- matrix(0, nrow=100, ncol=10)
for(j in 1:10) {
  set.seed (j)
  X [, j] <- rnorm(sampleSize)
}

# Generate response of binary problem with sum(X) > 0 -> 1 and 0 elsewhere
# with Gaussian noise
set.seed(-1)
error <- rnorm (100)
y <- ifelse((rowSums(X) + error) > 0, 1, 0)

# Calculate loss function with parameters (D=10, sigma=1, lambda=0)
# in one layer
calcLoss <- lossKDSN (parOpt=c(10, 1, 0), y=y, X=X, 
gammaPar=1, seedW=0)
str(calcLoss)

# Calculate loss function with parameters 
# (D=10, sigma=1, lambda=0.1, D=5, sigma=2, lambda=0.01) in two layers
calcLoss <- lossKDSN (parOpt=c(10, 1, 0.1, 5, 2, 0.01), 
y=y, X=X, gammaPar=1, seedW=rep(0, 2))
str(calcLoss)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line