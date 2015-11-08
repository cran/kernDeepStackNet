\name{fineTuneKDSN}
\alias{fineTuneKDSN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fine tuning of random weights of a given KDSN model
}
\description{
Weight matrices are randomly drawn given a KDSN model structure with prespecified number of levels, dimensions of the random Fourier transformations, precision parameter sigma and regularization parameter lambda. This function only supports the generalized cross-validation score. The model with the lowest gcv score is chosen.
}
\usage{
fineTuneKDSN(estKDSN, y, X, gammaPar=1, fineTuneIt=100, info=TRUE, seedInit=NULL, ncpus=1)
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
\item{gammaPar}{
Weight parameter in the generalized cross-validation score. Default is 1 (numeric scalar).
}
\item{fineTuneIt}{
Number of random generated iterations. Default is 100 (numeric scalar).
}
\item{info}{
Should additional informations be displayed? Default is TRUE (logical scalar).
}
\item{seedInit}{
Gives the seed initialization of the random seeds for weight matrix generation. Each element of the integer vector corresponds to one level. Default is NULL \code{\link[base]{Random}}.
}
\item{ncpus}{
Gives the number of cpus (threads), which should be used. The parallization is based on the parallel package. 
}
}
%\details{
%}
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
This procedure is already included in the other tuning algorithms. For further improving predictive accuracy, this function can be used after finding the structure of the KDSN. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{tuneKDSN}}, \code{\link{tuneLevelKDSN}}, \code{\link{tuneMboKDSN}}, \code{\link{tuneMboLevelKDSN}}
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
             sigma=c(0.5, 1, 2), lambda=c(1, 0.1, 0.01), 
             alpha=rep(0, 3), info=TRUE)

# Apply additional fine tuning
fitKDSN2 <- fineTuneKDSN(estKDSN=fitKDSN1, y=y, X=X, gammaPar=1, 
fineTuneIt=100, info=TRUE)

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