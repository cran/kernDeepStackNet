\name{predLogProb}
\alias{predLogProb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predictive logarithmic probability of Kriging model
}
\description{
Calculates the predictive logarithmic probability of Kriging model given the model parameters (likelihood function). It is used to choose better meta models in model-based optimization using Kriging. The observation used for prediction is left out in the training data and then estimated. 
}
\usage{
predLogProb(predMean, predSigma, y, X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{predMean}{
Predicted leave one out mean of the Kriging model of all observations (numeric vector).
}
\item{predSigma}{
Predicted leave one out variance of the Kriging model of all observations.
}
  \item{y}{
Numeric response vector of the outcome. Should be formated as a one column matrix.
}
  \item{X}{
Numeric design matrix of the covariates. All factors have to be prior encoded.
}
}
\details{
Gaussian processes are conditionaly normal distributed and therefore the normal likelihood is used. The leave out parameters can be efficiently computed without performing leave one out cross-validation.
}
\value{
Numeric Value of log-likelihood with leave-one-out parameters.
}
\references{
Carl Edward Rasmussen and Christopher K. I. Williams, (2006), 
\emph{Gaussian Processes for Machine Learning}
Massachusetts Institute of Technology
}
\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}
}

\note{
This function is not intended to be called directly by the user. Should only be used by experienced users, who want to customize the model. It is called in the model selection process of the kernel deep stacking network, e.g. \code{\link{tuneMboKDSN}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%}

\examples{
library(DiceKriging)
# Generate design of experiments
design.fact <- expand.grid(x1=seq(0,1,length=4), x2=seq(0,1,length=4))
y <- apply(design.fact, 1, branin) 

# Estimate Kriging model
km_fit <- km(design=design.fact, response=y, control=list(trace=FALSE), 
nugget.estim=TRUE, iso=FALSE)

# Calculate leave one out parameters and performance measure
loo <- leaveOneOut.km(model=km_fit, type="UK", trend.reestim=TRUE)
predLogProbs <- predLogProb(predMean=loo$mean, predSigma=loo$sd^2, y=y, X=X)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line