\name{tuneMboSharedSubsetKDSN}
\alias{tuneMboSharedSubsetKDSN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tuning subsets of KDSN with efficient global optimization and shared hyperparameters (experimental)
}
\description{
Estimates ensembles of (sparse) KDSNs with shared tuning parameters (still experimental). The number of levels is incorporated into the tuning process. Each part of the ensemble is tuned using model-based optimization as implemented in function \code{\link{tuneMboSharedCvKDSN}}. 
}
\usage{
tuneMboSharedSubsetKDSN (noSubsets=2, noSubsetRep=1, subSeed=NULL, y, X, alphaShared=1,
                         nStepMult=20, designMult=10,
                         lossFunc=devStandard, GenSAmaxCall=100,
                         varSelectShared=TRUE, dropHiddenShared=TRUE,
                         standX=TRUE, standY=FALSE, timeAlloc="constant", 
                         varPreSelect=TRUE, varPreSelpopSize=100, 
                         varPreSelMaxiter=100, maxLevels=10, nCores=1,
                         addInfo=1, saveOnDisk=FALSE, 
                         dirNameDisk=paste(tempdir(), "/ensembleModel", sep=""),
                         useAllSub=TRUE, trainPartition=0.5, noPartion=1,
                         EItype="EQI")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{noSubsets}{
Number of training data splits (integer scalar). Training data will be randomly split into disjoint subsets. Each subset will be used as training set and the remaining data as test set.
}
\item{noSubsetRep}{
Number of independent random subset draws (integer scalar).
}
\item{subSeed}{
Random seed of subset generation \code{\link[base]{Random}}. Default is NULL.
}
\item{y}{
Response matrix with one column.
}
\item{X}{
Design matrix. All factors must be already encoded.
}
\item{alphaShared}{
Weight parameter between lasso and ridge penalty (numeric vector) of each level. Default=0 corresponds to ridge penalty and 1 equals lasso.
}
\item{nStepMult}{
Multiplier, which affects how many steps the EGO algorithm is run, depending on the number of parameters to estimate.
}
\item{designMult}{
Multiplier, which affects how many initial design points are evaluated in the loss function, depending on the number of parameters to estimate.
}
\item{lossFunc}{
Specifies how the loss on the test data should be evaluated. Defaults to predictive deviance \code{\link{devStandard}}.
}
\item{GenSAmaxCall}{
Maximum number of function calls per parameter to estimate in generalized, simulated annealing. Higher values result in more accurate estimates, but the optimization process is slowed.
}
\item{varSelectShared}{
Specifies, if variables should be preselected by using randomized, dependence coefficient. Default is no variable selection. This setting is shared across all levels.
}
\item{dropHiddenShared}{
Should dropout be applied on the random Fourier transformation? Each entry corresponds to the one level. Default is without dropout (logical vector).
}
\item{standX}{
Should the design matrix be standardized by median and median absolute deviation? Default is TRUE.
}
\item{standY}{
Should the response be standardized by median and median absolute deviation? Default is FALSE.
}
\item{timeAlloc}{
Specifies how the new noise variance is influenced by iteration progress. Default is to use "constant"" allocation. The other available option is to specify "zero", which sets the future noise variance always to zero.
}
\item{varPreSelect}{
Should variables be pre-selected using RDC and genetic algorithm? Default is no. May consume a lot of start up time.
}
\item{varPreSelpopSize}{
Population size of the genetic algorithm (integer scalar).
}
\item{varPreSelMaxiter}{
Maximum number of generations of the genetic algorithm (integer scalar).
}
\item{maxLevels}{
Maximum number of levels possible to tune. Lower number speeds up tuning, but is less flexible. Default is to use 10 levels. (integer scalar)
}
\item{nCores}{
Number of threads to use in implizit calculation based on \code{\link[parallel]{parallel-package}} package. Default is serial processing. (integer scalar)
}
\item{addInfo}{
Should the progress during tuning be printed? (integer scalar) Value zero means no additional printing. The next step is addInfo=1 overall progress will be printed. addInfo=2 additionally prints progress of inner tuning procedures. Default is printing of overall progress. 
}
\item{saveOnDisk}{
Should the estimated models be saved on disk instead of memory? (logical scalar). If the data is high dimensional, the workspace memory may not be sufficient to store all of them at once. 
}
\item{dirNameDisk}{
Gives the directory and file name of the tuned SKDSN models. The Number at the end of the filename represents the subset. Default is the temporary R-directory.
}
\item{useAllSub}{
The data is split into \code{noSubsets} equal parts. If \code{useAllSub==TRUE} the tuning process is applied to all training parts. Otherwise the complete data set is split into one training and one validation set. The proportion is given by the argument \code{trainPartition}.
}
\item{trainPartition}{
Gives the proportion of the complete data set used in training. Only applied, if \code{useAllSub==FALSE}.
}
\item{noPartion}{
Gives the number of randomly drawn partitions of the complete data set. Only applied, if \code{useAllSub==FALSE}.
}
\item{EItype}{
Defines the type of the improvement criterion. The default \code{EQI} corresponds to the expected quantile improvement. As an alternative \code{EI} expected improvement is also possible.
}
}
\details{
The tuning parameters fixed across all levels. This leads to more parsimonious models and faster tuning. For additional flexibility the number of levels is not given in advance and also considered in tuning. 
}

\value{
Gives the best tuned kernel deep stacking network of class \code{KDSNensemble} or \code{KDSNensembleDisk}. For further details see \code{\link{predict.KDSNensemble}}, \code{\link{predict.KDSNensembleDisk}} and \code{\link{fitKDSN}}.
}

\references{
David Lopez-Paz and Philipp Hennig and Bernhard Schoelkopf, (2013), 
\emph{The Randomized Dependence Coefficient},
Max Planck Institute for Intelligent Systems, Germany

Victor Picheny, David Ginsbourger, Yann Richet, (2012),
\emph{Quantile-based optimization of Noisy Computer Experiments with Tunable Precision},
HAL-archives-ouvertes.fr, hal-00578550v3
}

\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}
}

%%\note{
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fitKDSN}}, \code{\link{tuneMboSharedCvKDSN}}, \code{\link[DiceKriging]{km}}, \code{\link[DiceKriging]{leaveOneOut.km}}, \code{\link[lhs]{maximinLHS}}, \code{\link{mboAll}}, \code{\link{mbo1d}}
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

# Draw cv training indices
library(caret)
cvTrainInd <- createFolds(y=y, k = 2, list = TRUE, returnTrain = TRUE)

# Define loss function
defLossFunc <- function(preds, ytest) {-c(auc(response=ytest, predictor=c(preds)))}

# Tune kernel deep stacking network by auc on test data
\dontrun{
tuned_KDSN_EGO_level <- tuneMboSharedCvKDSN (y=y, X=X, 
nStepMult=2, designMult=3, lossFunc=defLossFunc,
GenSAmaxCall=10, varPreSelpopSize=10, varPreSelMaxiter=10)
preds <- predict(tuned_KDSN_EGO_level, newx=Xtest)
library(pROC)
auc(response=ytest, predictor=c(preds))
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
