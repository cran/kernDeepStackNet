\name{tuneMboLevelGcvKDSN}
\alias{tuneMboLevelGcvKDSN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tuning of KDSN with efficient global optimization given level by cross-validation
}
\description{
Implements the efficient global optimization algorithm based on Kriging for kernel deep stacking networks (KDSN). This function uses generalized cross validation approximation of out of sample error without using external data. Preselection of variables is based on the randomized dependence coefficient (RDC). 
}
\usage{
tuneMboLevelGcvKDSN(y, X, levels=1, alpha=rep(0, levels), fineTuneIt=100, 
                    nStepMult=20, designMult=10, 
                    dimMax=round(sqrt(dim(X)[1])/2), addInfo=TRUE,
                    maxRunsMult=1, repMult=1, 
                    tol_input=.Machine$double.eps^0.25, 
                    EIopt="1Dmulti", GenSAmaxCall=100,
                    varSelect=rep(FALSE, levels), rdcRep=1, 
                    dropHidden=rep(FALSE, levels),
                    standX=TRUE, standY=FALSE, timeAlloc="constant", 
                    varPreSelect=FALSE, varPreSelpopSize=100, varPreSelMaxiter=100,
                    EItype="EQI")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{y}{
Response matrix with one column.
}
\item{X}{
Design matrix. All factors must be already encoded.
}
\item{levels}{
Maximum number of levels for the kernel deep stacking network (integer scalar).
}
\item{alpha}{
Weight parameter between lasso and ridge penalty (numeric vector) of each level. Default=0 corresponds to ridge penalty and 1 equals lasso.
}
\item{fineTuneIt}{
Number of drawn random weight matrices in fine tuning (integer scalar). If set to zero, no fine tuning is done.
}
\item{nStepMult}{
Multiplier, which affects how many steps the EGO algorithm is run, depending on the number of parameters to estimate.
}
\item{designMult}{
Multiplier, which affects how many initial design points are evaluated in the loss function, depending on the number of parameters to estimate.
}
\item{dimMax}{
Maximal dimension of the random Fourier transformation. The effective number of parameters is dimMax*2. The default heuristic depends on the sample size.
}
\item{addInfo}{
Should additional information be printed during estimation? Default is TRUE.
}
\item{maxRunsMult}{
Multiplies the base number of iterations in the conditional one dimensional optimization. Default is to use the number of hyperparameters. See \code{\link{optimize1dMulti}}.
}
\item{repMult}{
Multiplies the base number of random starting values in the conditional one dimensional optimization to avoid local optima. Default is to use the number of hyperparameters. See \code{\link{optimize1dMulti}}.
}
\item{tol_input}{
Convergence criteria of each one dimensional sub-optimization. Higher values will be more accurate, but require much more function evaluations. Default is the fourth root of the machine double accuracy. See \code{\link{optimize1dMulti}}.
}
\item{EIopt}{
Specifies which algorithm is used to optimize the expected improvement criterion. Two alternatives are available "1Dmulti" and "GenSA". The former uses the conditional 1D algorithm and the latter generalized, simulated annealing.
}
\item{GenSAmaxCall}{
Maximum number of function calls per parameter to estimate in generalized, simulated annealing. Higher values result in more accurate estimates, but the optimization process is slowed.
}
\item{varSelect}{
Specifies, if variables should be preselected by using randomized, dependence coefficient. Default is no variable selection in all levels.
}
\item{rdcRep}{
Number of repetitions for the randomized dependence coefficient \code{\link{rdcVarOrder}}.
}
\item{dropHidden}{
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
\item{EItype}{
Defines the type of the improvement criterion. The default \code{EQI} corresponds to the expected quantile improvement. As an alternative \code{EI} expected improvement is also possible.
}
}
%\details{
%The direct tuning algorithm (see for reference \code{\link{tuneLevelKDSN}}) is adapted %to surrogate model optimization. The algorithm is similiar and uses a model-based %optimization procedure.
%}

%return(list(par=c(ego_result), value=attr(ego_result, "funcVal")))

\value{
Gives the best tuned kernel deep stacking network of class \code{k_DSN_rft} given a specific level (see \code{fitKDSN}).
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

\note{
This function is supplied to optimize units KDSN with many levels. Complex data generating mechanisms may be more efficiently represented with large number of levels. The computation time increases progressive the more levels are added. Higher values of tol_input, fineTuneIt, maxRuns, repetitions may increase performance. The computation time of the tuning algorithm increases mostly with higher values of dimMax. \\
\\
The variable selection \code{varSelect=TRUE} calculates the rdc of all variables. Then the variables above a pre-specified threshold are choosen as important and all others are discarded. In the estimation process, in each level the least important variable according to rdc, is removed from the design matrix. Only observed variables may be removed. Latent variables are always used in the model.  
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[DiceKriging]{km}}, \code{\link[DiceKriging]{leaveOneOut.km}}, \code{\link[lhs]{maximinLHS}}, \code{\link{mboAll}}, \code{\link{mbo1d}}
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

# Tune kernel deep stacking network by auc on test data
\dontrun{
tuned_KDSN_EGO_level <- tuneMboLevelGcvKDSN (y=y, X=X, 
levels=2, fineTuneIt=10, nStepMult=2, designMult=5)
preds <- predict(tuned_KDSN_EGO_level, newx=Xtest)
library(pROC)
auc(response=ytest, predictor=c(preds))
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line