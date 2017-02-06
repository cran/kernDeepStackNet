\name{lossApprox}
\alias{lossApprox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kernel deep stacking network loss function
}
\description{
Computes the generalized cross-validation (gcv) score of a tuning parameter configuration given the data. The tuning parameters consists of dimension of the random Fourier transform (=D), the variance of the Gaussian kernels (=sigma) and the regularization parameter of the kernel ridge regression (=lambda). It is a wrapper around the function \code{\link{lossGCV}} and therefore only supports the generalized cross-validation score. First a kernel deep stacking network is fitted. Then the inputs are prepared to calculate the gcv score. 
}
\usage{
lossApprox(parOpt, y, X, levels, seedW=NULL,
           varSelect=rep(FALSE, levels), varRanking=NULL, 
           alpha=rep(0, levels), dropHidden=rep(FALSE, levels),
           seedDrop=NULL, standX=TRUE, standY=FALSE,
           namesParOpt=rep(c("dim", "sigma", "lambdaRel"), levels), 
           gammaPar=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parOpt}{
The numeric parameter vector of the tuning parameters. The length of the vector is required to be divisible with three without rest. The order of the parameters is (D, sigma, lambda) in one level. Each additional level is added as new elements, e.g. with two levels (D_1, sigma_1, lambda_1, D_2, sigma_2, lambda_2), etc.
}
  \item{y}{
Numeric response vector of the outcome. Should be formated as a one column matrix.
}
  \item{X}{
Numeric design matrix of the covariates. All factors have to be prior encoded.
}
\item{levels}{
  Number of levels of the kernel deep stacking network (integer scalar).
}
  \item{seedW}{
Seeds for the random Fourier transformation (integer vector). Each value corresponds to the seed of one level. The length must equal the number of levels of the kernel deep stacking network.
}
\item{varSelect}{
Should unimportant variables be excluded in each level? Default is that all available variables are used. (logical scalar)
}
\item{varRanking}{
Defines a variable ranking in increasing order. The first variable is least important and the last is the most important. (integer vector)
}
\item{alpha}{
Weight parameter between lasso and ridge penalty (numeric vector) of each level. Default=0 corresponds to ridge penalty and 1 equals lasso.
}
\item{dropHidden}{
Should dropout be applied on the random Fourier transformation? Each entry corresponds to the one level. Default is without dropout (logical vector).
}
\item{seedDrop}{
Specifies the seed of the random dropouts in the calculation of random Fourier transformation per level. Default is random (integer vector).
}
\item{standX}{
Should the design matrix be standardized by median and median absolute deviation? Default is TRUE.
}
\item{standY}{
Should the response be standardized by median and median absolute deviation? Default is FALSE.
}
\item{namesParOpt}{
Gives the names of the argument parOpt (character vector). It is used to encode the parameters into the correct structure suitable in fitting the kernel deep stacking network.
}
  \item{gammaPar}{
Weighting parameter (numeric scalar), which specifies how the generalized cross-validation score is penalized by the effective degrees of freedom. Default value is 1.
}
}
%\details{
%
%}
\value{
Numeric scalar with gives the generalized cross-validation score. The kernel deep stacking network used to calculate the score is available in the attributes.
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
This function is not intended to be called directly by the user. Should only be used by experienced users, who want to customize the model. %It is called in the model selection process of the kernel deep stacking network, e.g. \code{\link{tuneKDSN}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lossCvKDSN}}, \code{\link{lossGCV}}, \code{\link{lossSharedCvKDSN}}, 
\code{\link{lossSharedTestKDSN}}
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
calcLoss <- lossApprox (parOpt=c(10, 1, 0), y=y, X=X, 
levels=1, seedW=0)
str(calcLoss)

# Calculate loss function with parameters 
# (D=10, sigma=1, lambda=0.1, D=5, sigma=2, lambda=0.01) in two layers
calcLoss <- lossApprox (parOpt=c(10, 1, 0.1, 5, 2, 0.01), 
y=y, X=X, levels=1, seedW=0)
str(calcLoss)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line