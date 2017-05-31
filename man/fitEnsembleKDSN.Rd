\name{fitEnsembleKDSN}
\alias{fitEnsembleKDSN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit an ensemble of KDSN (experimental)
}
\description{
Given a KDSN structure, either tuned \code{\link{tuneMboLevelCvKDSN}}, \code{\link{tuneMboSharedCvKDSN}} or estimated \code{\link{fitKDSN}}, bootstrap samples are choosen of the original data. Then new random Fourier weights and random Dropout (if specified) are generated. A prespecified number of ensembles is fitted.
}
\usage{
fitEnsembleKDSN (estKDSN, y, X, ensembleSize=100, 
                 bagging=rep(FALSE, ensembleSize), seedBag=NULL, 
                 randSubsets=rep(FALSE, ensembleSize),
                 seedRandSubSize=NULL, seedRandSubDraw=NULL,
                 seedW1=NULL, seedW2=NULL, 
                 seedDrop1=NULL, seedDrop2=NULL, 
                 info=TRUE, nCores=1, saveOnDisk=FALSE, 
                 dirNameDisk=paste(tempdir(), "/ensembleModel", sep=""))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{estKDSN}{
Model of class "KDSN" \code{\link{fitKDSN}}.
}
\item{y}{
Response of the regression (must be in one column matrix format).
}
\item{X}{
Design matrix. All factors must be already encoded.
}
\item{ensembleSize}{
Number of ensembles. Default is 100 (integer scalar).
}
\item{bagging}{
Should bagging be used in ensemble? Default is FALSE. Each entry corresponds to usage in one ensemble (logical vector).
}
\item{seedBag}{
Gives the seed initialization of the bootstrap samples. Each element of the integer vector corresponds to one ensemble. Default is NULL \code{\link[base]{Random}}.
}
\item{randSubsets}{
Should random subsets be used in ensemble? Default is FALSE. Each entry corresponds to usage in one ensemble (logical vector).
}
\item{seedRandSubSize}{
Seed for size of random subsets selection. Each entry corresponds to the seed of one ensemble (integer vector).
}
\item{seedRandSubDraw}{
Seed for random subsets selection given size. Each entry corresponds to the seed of one ensemble (integer vector).
}
\item{seedW1}{
Gives the seed initialization of the bootstrap samples. The random number generation is split into two parts: This seed represents drawing random integer numbers. Each element of the integer vector corresponds to one ensemble. Default is NULL \code{\link[base]{Random}}.
}
\item{seedW2}{
Gives the seed initialization of the bootstrap samples. The random number generation is split into two parts: This seed represents drawing random signs. Each element of the integer vector corresponds to one ensemble. Default is NULL \code{\link[base]{Random}}.
}
\item{seedDrop1}{
Gives the seed initialization of the dropout applied in the random Fourier transformation. The random number generation is split into two parts: This seed represents drawing random integer numbers. Each element of the integer vector corresponds to one ensemble. Default is NULL \code{\link[base]{Random}}.
}
\item{seedDrop2}{
Gives the seed initialization of the dropout applied in the random Fourier transformation. The random number generation is split into to parts: This seed represents drawing random signs. Each element of the integer vector corresponds to one ensemble. Default is NULL \code{\link[base]{Random}}.
}
\item{info}{
Should additional informations be displayed? Default is TRUE (logical scalar).
}
\item{nCores}{
Gives the number of threads, which will be executed by the \code{\link[parallel]{parallel-package}} package (integer scalar). Defaults to no parallel processing.
}
\item{saveOnDisk}{
Should all models be saved on hard disk instead of in the workspace? (Logical scalar) If the data sets are big and the KDSN is deep, then one model will occupy lots of workspace. Sometimes it is prohibitive to store the complete ensemble into workspace. In this approach every ensemble model is saved in an temporary Folder and the names of the files are saved. Later in prediction only one model is stored in workspace and stored until predictions are made. Only supported for single thread execution. 
}
\item{dirNameDisk}{
Only relevant if saveOnDisk=TRUE. Gives the path, where the fitted models will be saved.
}
}
\details{
Note that this function is still experimental. It is important that the estimated KDSN has a custom set seed (argument "estKDSN") to ensure reproducibility. Otherwise the model is refitted using the default seed values 
seq(0, (Number of Levels-1), 1) are used for random Fourier transformation. 

The disadvantage on using argument saveOnDisk=TRUE is, that saving and compressing the objects takes more time instead of using them together into workspace.
}
\value{
Model of class "KDSNensemble" or "KDSNensembleDisk" if saveOnDisk=TRUE. 
For reproduceability all seeds of the random iterations are stored as attributes.
}
%\references{
%Simon N. Wood, (2006), 
%\emph{Generalized Additive Models: An Introduction with R},
%Taylor \& Francis Group LLC
%}
\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}
}

%\note{
%This procedure is already included in the other tuning algorithms. For further %improving predictive accuracy, this function can be used after finding the structure of %the KDSN. 
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%\code{\link{tuneMboCvKDSN}}, \code{\link{tuneMboLevelCvKDSN}}
%}

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
# Initial seed should be supplied in fitted model!
fitKDSN1 <- fitKDSN(y=y, X=X, levels=3, Dim=c(20, 10, 5), 
             sigma=c(0.5, 1, 2), lambdaRel=c(1, 0.1, 0.01), 
             alpha=rep(0, 3), info=TRUE, seedW=c(0, 1:2))

# Fit 10 ensembles
fitKDSNens <- fitEnsembleKDSN (estKDSN=fitKDSN1, y=y, X=X, ensembleSize=10, 
                 seedBag=1:10, seedW1=101:110, seedW2=-(101:110), 
                 seedDrop1=3489:(3489+9), seedDrop2=-(3489:(3489+9)), info=TRUE)

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
auc1 <- auc(response=ytest, predictor=c(preds))
predsMat <- predict(fitKDSNens, Xtest)
preds <- rowMeans(predsMat)
auc2 <- auc(response=ytest, predictor=c(preds))
auc1 < auc2 # TRUE
# The ensemble predictions give a better test auc than the original model
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line