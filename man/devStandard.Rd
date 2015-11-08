\name{devStandard}
\alias{devStandard}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predictive deviance of a linear model
}
\description{
Calculates the deviance on new data observations. The predictive deviance measures how far the predicted values are apart from the saturated model.
}
\usage{
devStandard(preds, ytest)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{preds}{
Predictions of the specified model (numeric vector).
}
\item{ytest}{
Data values of the response in the test data.
}
}
\details{
In the "Gaussian" case it is defined to be the residual sum of squares. \code{ytest} are the test observations and \code{preds} are the predicted values of the model on the test data.
}
\value{
Predictive deviance of the linear model, given predictions and test data (numeric scalar).
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
This function is not intended to be called directly by the user. Should only be used by experienced users, who want to customize the model. It is called in the model selection process of the kernel deep stacking network with cross-validation, e.g. \code{\link{lossCvKDSN}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%}

\examples{
##################
# Fit Gaussian glm

set.seed(10)
x <- matrix(rnorm(100*20),100,20)
set.seed(100)
y <- rnorm(100)
fit1 <- glm(formula=y ~ ., data=data.frame(x))
preds <- predict(fit1, type="response")
# Performance on training data
all.equal(devStandard(preds=preds, ytest=y), fit1$deviance)
# Performance on random test data
set.seed(100)
yTest <- simulate(fit1)
devStandard(preds=preds, ytest=yTest)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line