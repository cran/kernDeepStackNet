\name{tuneLevelGridKDSN}
\alias{tuneLevelGridKDSN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tune kernel deep stacking network by direct optimization over a set of levels
}
\description{
Supplies direct optimization of the tuning parameters of the kernel deep stacking network for a given level. Direct optimization means, that the generalized cross-validation score is minimized directly by the one dimensional optimization procedure. The optimization procedure is conditional on the set of levels. For each level the direct optimization regarding the GCV score is performed. The tuned KDSN with the smallest GCV score is returned.
}
\usage{
tuneLevelGridKDSN(y, X, maxRuns=3, repetitions=5, levelSet, gammaPar=1, 
                  fineTuneIt=100, tol_input=.Machine$double.eps^0.25, addInfo=TRUE, 
                  dimMax=round(sqrt(dim(X)[1])/2), nCoresInner=1, nCoresOuter=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{y}{
Numeric matrix of the response vector with one column.
}
\item{X}{
Design matrix of the covariates. All factors must be already encoded. Does not need an intercept column.
}
\item{maxRuns}{
Maximum number of iterations of the one dimensional algorithm.
}
\item{repetitions}{
Maximum number of random starting values of the one dimensional algorithm.
}
\item{levelSet}{
Gives a set of levels to perform grid search (integer vector).
}
\item{gammaPar}{
Weighting parameter in the generalized cross-validation score.
}
\item{fineTuneIt}{
Maximum number of fine tuning operations. Fine tuning consists of randomly drawn weights of the random Fourier transformation (integer scalar).
}
\item{tol_input}{
Accuracy threshold of the one dimensional algorithm (numeric scalar).
}
\item{addInfo}{
Should additional information be displayed during optimization? (logical value). Default is TRUE.
}
\item{dimMax}{
Maximal dimension of the random Fourier transformation. The effective number of parameters is dimMax*2. The default heuristic depends on the sample size.
}
\item{nCoresInner}{
Number of threads for computing initial designs of direct and fine tuning in parallel mod. The cluster ist set up using package parallel. Default is serial mode with nCoresInner=1.
}
\item{nCoresOuter}{
Number of threads for computing each MBO tuning in parallel mode. The cluster ist set up using package parallel. Default is serial mode with nCoresOuter=1.
}
}
\details{
First the starting values and bounds of all tuning parameters are set: 
\itemize{
\item{Dim: } {Bounds are set from 1 to number of observations n. Starting value is \code{round((dimMax+1)/2)}
\item{sigma: } {The maximum and minimum of the squared distances of the scaled design matrix are used as bounds. The 50 \% Quantile is the starting value}
\item{lambda: } {Bounds are the interval [0, 10]. Start is set to 1}
}
Then the count of levels is set to 1 and fix the random Fourier transformation weights to a specific seed. Given one level, the algorithm finds the best tuning parameters of the kernel deep stacking network. After this a fine tuning of the random weights is done, by randomly drawing network configurations. The best network is choosen according to the gcv score.

Repeat the same with one additional level and set the best previous tuning parameters as starting values in the lower level. Compare the best solutions of both level configurations and stop if the next level yields no additional benefit. Otherwise continue, until there is no improvement or the maximum number of levels is reached.
}
}

\value{
Gives the best tuned kernel deep stacking network of class \code{KDSN} given a specific level.
}

\section{Warning }{
Do not specify nCoresInner>1 and nCoresOuter>1, because the threads would conflict. Either choose the inner or outer loops to be parallelised. 
}

%\references{
%}

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
This function is supplied to optimize units KDSN with many levels. Very complex data generating mechanisms may be more efficiently represented with large number of levels. The computation time increases progressive the more levels are added. Higher values of tol_input, fineTuneIt, maxRuns, repetitions will increase performance considerably. The computation time of the tuning algorithm increases mostly with higher values of dimMax.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fitKDSN}}, \code{\link{optimize1dMulti}}, \code{\link{lossKDSN}}, \code{\link{randomFourierTrans}}, \code{\link{tuneKDSN}}, \code{\link{tuneLevelKDSN}}
}

\examples{
# Generate small sample of 20 observations of a binary classification task
# Due to keeping the example as fast as possible, the parameters of the tuning 
# algorithm are set for low accuracy. Higher values of tol_input, fineTuneIt, 
# maxRuns, repetitions will increase performance considerably.
library(pROC)

# Generate design matrix
sampleSize <- 20
X <- matrix(0, nrow=sampleSize, ncol=5)
for(j in 1:5) {
  set.seed (j)
  X [, j] <- rnorm(sampleSize)
}

# Generate response of binary problem with sum(X) > 0 -> 1 and 0 elsewhere
set.seed(-1)
error <- rnorm (sampleSize)
y <- ifelse((rowSums(X) + error) > 0, 1, 0)

# Generate test data
Xtest <- matrix(, nrow=sampleSize, ncol=5)
for(j in 1:5) {
  set.seed (j*2+1)
  Xtest [, j] <- rnorm(sampleSize)
}

# Generate test response
set.seed(-10)
error <- rnorm (sampleSize)
ytest <- ifelse((rowSums(Xtest) + error) > 0, 1, 0)

# Tune kernel deep stacking network
\dontrun{
tuned_KDSN_level <- tuneLevelGridKDSN (y=y, X=X, 
maxRuns=2, repetitions=2, levelSet=c(1, 2), gammaPar=1, fineTuneIt=10, 
tol_input=.Machine$double.eps^0.05)
preds <- predict(tuned_KDSN_level, newx=Xtest)
library(pROC)
auc(response=ytest, predictor=c(preds))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line