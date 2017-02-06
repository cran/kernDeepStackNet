\name{fitKDSN}
\alias{fitKDSN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit kernel deep stacking network with random Fourier transformations
}
\description{
Estimates the kernel deep stacking network with identity link function (Gaussian) for the response in each level. Sparsity may be included with the addition of dropout, variable pre- and internal selection as well as L1 penalty of the cofficients. 
}
\usage{
fitKDSN(y, X, levels=1, Dim=rep(round(sqrt(dim(X) [1]) / 2), levels), 
             sigma=rep(1, levels), lambdaRel=rep(0.5, levels), 
             alpha=rep(0, levels), info=FALSE, seedW=NULL, standX=TRUE, standY=FALSE,
             varSelect=rep(FALSE, levels), varRanking=NULL, varExFreq=rep(NA, levels),
             dropHidden=rep(FALSE, levels), 
             dropHiddenProb=rep(NA, levels), 
             seedDrop=NULL, baggingInd=NULL, randSubsetInd=NULL,
             maxItOpt=10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{y}{
Response variable as numeric vector.
}
\item{X}{
Design matrix. All factors must be encoded, e.g. dummy coding.
}
\item{levels}{
Number of levels of the kernel deep stacking network (integer scalar).
}
\item{Dim}{
Dimension of the random Fourier transformations in each level (integer vector). The first entry corresponds to the first level, the second entry to the second level and so on.
}
\item{sigma}{
Variance of the Gaussian kernel of each level. The higher the variance, on average the more different transformed observations will be generated.
}
\item{lambdaRel}{
Relative ridge regularization parameter of each level. The relative ridge regularization parameter is limited to the interval [0, 1] and 0 corresponds to no shrinkage and 1 stands for maximum penalty to shrink all coefficients except the intercept to zero. Default is set in the middle between those extremes. 
}
\item{alpha}{
Weight parameter between lasso and ridge penalty (numeric vector) of each level. Default=0 corresponds to ridge penalty and 1 equals lasso.
}
\item{info}{
Determines if additional infos of the level computation are displayed during training (logical value). Default is FALSE.
}
\item{seedW}{
Random seed for drawing Fourier transformation weights of the multivariate normal distribution (integer vector). Each entry in the vector corresponds to one level.
}
\item{standX}{
Should the design matrix be standardized by median and median absolute deviation? Default is TRUE.
}
\item{standY}{
Should the response be standardized by median and median absolute deviation? Default is FALSE.
}
\item{varSelect}{
Should unimportant variables be excluded in each level? Default is that all available variables are used. (logical scalar)
}
\item{varRanking}{
Defines a variable ranking in increasing order. The first variable is least important and the last is the most important. (integer vector)
}
\item{varExFreq}{
Gives the relative frequency of variables to drop within the interval [0, 1]. Latent variables are always included in the model.
}
\item{dropHidden}{
Should dropout be applied on the random Fourier transformation? Each entry corresponds to the one level. Default is without dropout in all levels (logical vector).
}
\item{dropHiddenProb}{
Probability of inclusion of cells in random Fourier transformation. Each entry corresponds to the probability of one level.
}
\item{seedDrop}{
Specifies the seed of the random dropouts in the calculation of random Fourier transformation per level. Default is random (integer vector).
}
\item{baggingInd}{
Gives the indices of the bootstrap samples in list format. Each entry of the list (integer vector) corresponds to one level. 
}
\item{randSubsetInd}{
Gives the indices of the random variable subsets in list format. Each entry of the list (integer vector) corresponds to one level. 
}
\item{maxItOpt}{
Gives the maximum number of iterations in the optimization procedure of the elastic net regression. For details see \code{\link[glmnet]{glmnet}}. 
}
}
\details{
A kernel deep stacking network is an approximation of artificial neural networks with multiple number of levels. In each level there is a random Fourier transformation of the data, based on Gaussian kernels, which can be represented as an infinite series. The random Fourier transform converges exponentially fast to the Gaussian kernel matrix by increasing dimension. Then a kernel ridge regression is applied given the transformed data. The predictions of the all last levels are included as covariates in the next level.
}
\value{
\itemize{
\item{Output: }{A list with objects:}
\itemize{
\item{rftX: } {List of random Fourier transformation matrices for each level}
\item{linWeights: } {List of estimated parameter vectors of the ridge regression for each level}
}
\item{Input: }{A list with all specified input parameter, e. g. number of levels, dimension of each level, ridge penalty, etc.}
}
}

\references{
Po-Seng Huang and Li Deng and Mark Hasegawa-Johnson and Xiaodong He, (2013), 
\emph{Random Features for kernel deep convex network},
Proceedings IEEE International Conference on Acoustics, Speech, and 
Signal Processing (ICASSP)

D.S. Broomhead and David Lowe, (1988), 
\emph{Radial Basis Functions, Multi-Variable Functional Interpolation and 
Adaptive Networks},
Controller HMSO, London

Jerome Friedman and Trevor Hastie and Rob Tibshirani, (2008), 
\emph{Regularization Paths for Generalized Linear Models via Coordinate Descent},
Department of Statistics, Stanford University
}

%References in article: Number 4, 20, 25, 15

%@book{RBFfuncApprox,
% title={Radial Basis Functions, Multi-Variable Functional Interpolation and Adaptive Networks},
% author={D.S. Broomhead and David Lowe},
% year={1988},
% publisher={Controller HMSO, London}
%}


\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}
}

\note{
It is recommended to standardize the design matrix before training. Otherwise the resulting values of the Fourier transform can be either all very small or very high, which results in poor prediction performance.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{predict.KDSN}}, \code{\link{randomFourierTrans}}, \code{\link{robustStandard}}, \code{\link[glmnet]{glmnet}} 
}

\examples{
####################################
# Example with binary outcome

# Generate covariate matrix
sampleSize <- 100
X <- matrix(0, nrow=100, ncol=5)
for(j in 1:5) {
  set.seed (j)
  X [, j] <- rnorm(sampleSize)
}

# Generate bernoulli response
rowSumsX <- rowSums(X)
logistVal <- exp(rowSumsX) / (1 + exp(rowSumsX))
set.seed(-1)
y <- sapply(1:100, function(x) rbinom(n=1, size=1, prob=logistVal[x]))

# Fit kernel deep stacking network with three levels
fitKDSN1 <- fitKDSN(y=y, X=X, levels=3, Dim=c(20, 10, 5), 
             sigma=c(6, 8, 10), lambdaRel=c(1, 0.1, 0.01), 
             alpha=rep(0, 3), info=TRUE, 
             seedW=c(1882335773, -1640543072, -931660653))

# Generate new test data
sampleSize <- 100
Xtest <- matrix(0, nrow=100, ncol=5)
for(j in 1:5) {
  set.seed (j+50)
  Xtest [, j] <- rnorm(sampleSize)
}
rowSumsXtest <- rowSums(Xtest)
logistVal <- exp(rowSumsXtest) / (1 + exp(rowSumsXtest))
set.seed(-1)
ytest <- sapply(1:100, function(x) rbinom(n=1, size=1, prob=logistVal[x]))

# Evaluate on test data with auc
library(pROC)
preds <- predict(fitKDSN1, Xtest)
auc(response=ytest, predictor=c(preds))

####################################
# Example with continuous outcome

# Generate covariate matrix
sampleSize <- 100
X <- matrix(0, nrow=100, ncol=5)
for(j in 1:5) {
  set.seed (j)
  X [, j] <- rnorm(sampleSize)
}

# Generate Gaussian random variable conditional on the covariates
linPred <- 1*X[,1] + 2*X[,2] + 0.5*X[,3] + 0.5
set.seed(-1)
y <- sapply(1:100, function(x) rnorm(n=1, mean=linPred[x], sd=2))

# Fit kernel deep stacking network with five levels
fitKDSN1 <- fitKDSN(y=y, X=X, levels=5, Dim=c(40, 20, 10, 7, 5), 
             sigma=c(0.125, 0.25, 0.5, 1, 2), lambdaRel=c(1, 0.5, 0.1, 0.01, 0.001), 
             alpha=rep(0, 5), info=TRUE, 
             seedW=c(-584973296, -251589341, -35931358, 313178052, -1322344272))

# Generate new data
sampleSize <- 100
Xtest <- matrix(0, nrow=100, ncol=5)
for(j in 1:5) {
  set.seed (j+10)
  Xtest [, j] <- rnorm(sampleSize)
}
linPred <- 1*Xtest[,1] + 2*Xtest[,2] + 0.5*Xtest[,3] + 0.5
set.seed(-10)
ytest <- sapply(1:100, function(x) rnorm(n=1, mean=linPred[x], sd=2))

# Predictions on new data and compute root mean squared error
preds <- predict(obj=fitKDSN1, newx=Xtest)
RMSE <- sqrt(mean((ytest-c(preds))^2))
RMSE
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line