\name{fineTuneCvKDSN}
\alias{fineTuneCvKDSN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fine tuning of random weights of a given KDSN model
}
\description{
Weight matrices are randomly drawn given a KDSN model structure with prespecified number of levels, dimensions of the random Fourier transformations, precision parameter sigma and regularization parameter lambda. This function supports arbitrary loss functions, which are evaluated on test data. The model with the lowest loss is chosen. 
}
\usage{
fineTuneCvKDSN(estKDSN, y, X, fineTuneIt=100, info=TRUE, 
seedInitW=NULL, seedInitD=NULL, ncpus=1, cvIndex, lossFunc=devStandard)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{estKDSN}{
Model of class "KDSN" \code{\link{fitKDSN}}
}
\item{y}{
Response of the regression (must be in one column matrix format).
}
\item{X}{
Design matrix. All factors must be already encoded.
}
\item{fineTuneIt}{
Number of random generated iterations. Default is 100 (numeric scalar).
}
\item{info}{
Should additional informations be displayed? Default is TRUE (logical scalar).
}
\item{seedInitW}{
Gives the seed initialization of the random seeds for weight matrix generation. Each element of the integer vector corresponds to one level. Default is NULL \code{\link[base]{Random}}.
}
\item{seedInitD}{
Gives the seed initialization of the random seeds for dropout generation. Each element of the integer vector corresponds to one level. Default is NULL \code{\link[base]{Random}}.
}
\item{ncpus}{
Gives the number of cpus (threads), which should be used. The parallization is based on the parallel package. 
}
\item{cvIndex}{
Index of cross-validation indices. The indices represent the training data. Must be supplied as list, the required format is identical to the output of the \code{\link[caret]{createFolds}} with argument returnTrain=TRUE. 
}
\item{lossFunc}{
Specifies how the loss on the test data should be evaluated. Defaults to predictive deviance \code{\link{devStandard}}.
}
}
\details{
It is important that the estimated KDSN has a custom set seed (argument "estKDSN") to ensure reproducibility. Otherwise the model is refitted using the default seed values 
seq(0, (Number of Levels-1), 1) are used for random Fourier transformation.
}
\value{
Model of class "KDSN". For reproduceability the best seeds of the random iterations are stored in the final model.
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
This procedure is already included in the other tuning algorithms, e. g. \code{\link{tuneMboLevelCvKDSN}}. It can be used after finding the structure of the KDSN. This may improve predictive accuracy in some cases, but it is not always necessary. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{tuneMboLevelCvKDSN}}, \code{\link{tuneMboSharedCvKDSN}}, \code{\link{tuneMboLevelGcvKDSN}}
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

# Generate test indices
library(caret)
set.seed(123)
cvFoldIndices <- createFolds(y=y, k=2, returnTrain=TRUE)

# Fit kernel deep stacking network with three levels
# Initial seed should be supplied in fitted model!
fitKDSN1 <- fitKDSN(y=y, X=X, levels=3, Dim=c(20, 10, 5), 
             sigma=c(0.5, 1, 2), lambdaRel=c(1, 0.1, 0.01), 
             alpha=rep(0, 3), info=TRUE, seedW=c(0, 1:2))

# Apply additional fine tuning based on predictive deviance
fitKDSN2 <- fineTuneCvKDSN(estKDSN=fitKDSN1, y=y, X=X, 
fineTuneIt=25, info=TRUE, cvIndex=cvFoldIndices)

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
preds <- predict(fitKDSN2, Xtest)
auc2 <- auc(response=ytest, predictor=c(preds))
auc1 < auc2 # TRUE
# The fine tuning improved the test auc
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line