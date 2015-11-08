\name{optimize1dMulti}
\alias{optimize1dMulti}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
One dimensional optimization of multivariate loss functions
}
\description{
Applies a stepwise one dimensional algorithm conditional on the previous best parameter values. The one dimensional algorithm is a combination of golden search and linear interpolation.
}
\usage{
optimize1dMulti(f_input, interval_matrix, maxRuns=3, repetitions=5, 
tol_input=.Machine$double.eps^0.25, x_0=NULL, addInfo=TRUE,
nCores=1, envir=parent.frame(), directUse=TRUE, OptTypePar="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{f_input}{
General loss function with multiple inputs and one numeric outcome value. The function must have a x as first argument, e.g. f(x).
}
\item{interval_matrix}{
Gives the bounds of the input parameters (numeric matrix). The dimension of the matrix is 2 x I and I is the number of inputs. The first row is defined as the minimum and the second row as the maximum bound of the parameter in the i-th column.
}
\item{maxRuns}{
Maximal number of iterations in the conditional optimization. Default is three.
}
\item{repetitions}{
Maximal number of random starting values to avoid local optima. Default is five.
}
\item{tol_input}{
Convergence criteria of each one dimensional sub-optimization. Higher values will be more accurate, but require much more function evaluations. Default is the fourth root of the machine double accuracy.
}
\item{x_0}{
Intitial parameter values. If not supplied, a random number is drawn between the given bounds.
}
\item{addInfo}{
Should the optimization print additional information? (logical value) Default is TRUE. The information consists of printing, which parameters have been optimized, which iteration and repetition the optimization process just finished.
}
\item{nCores}{
Specifies how many cores are used in tuning (integer scalar). Default=1 is serial processing.
}
\item{envir}{
Internal variable to store environments. Default is to look up the next higher level environment. Modification is unnecessary.
}
\item{directUse}{
Specifies if the optimization is carried out inside a higher level function (see ) or directly (logical scalar). Only relevant in parallel mode initialization.
}
\item{OptTypePar}{
Internal variable to store, which higher level function was called previously. Modification has no effect.
}
}
\details{
First the supplied initional value \code{x_0} is used. If not present a random intitialisation is used. Then each parameter is optimized, conditional on the other starting values. An iteration is finished, if all parameters have been optimized once. The procedure is repeated until convergence or \code{maxRuns} is reached. If \code{repetitions} is higher than 1, new starting values are generated and the process starts new. In the end the best parameter values according the lowest function value are returned.
}
\value{
\itemize{
\item{}{List with following components: }
\itemize{
\item{minimum: } {Best parameter values}
\item{objective: } {Minimum of function evaluated at the best parameter values}
}
}
}
%Output <- list(minimum=Results [[Index]]$minimum, objective=Results [[Index]]$objective)

\references{
R. Brent, (1973),
\emph{Algorithms for Minimization without Derivatives},
Englewood Cliffs N.J.: Prentice-Hall
}

%@book{optimizeAlg,
% title={Algorithms for Minimization without Derivatives},
% author={R. Brent},
% year={1973},
% publisher={Englewood Cliffs N.J.: Prentice-Hall}
%}

\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}
}

%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{tuneKDSN}}, \code{\link{tuneLevelKDSN}}, \code{\link{tuneMboKDSN}}, \code{\link{tuneMboLevelKDSN}}
}

\examples{
library(globalOptTests)
# Two dimensional task
getDefaultBounds("AluffiPentini")
oneDres <- optimize1dMulti (f_input=function (x) 
goTest (x, fnName="AluffiPentini", checkDim=TRUE), 
maxRuns=3, repetitions=5, interval_matrix=matrix(rep(c(-12, 10), 2), ncol=2), 
tol_input=.Machine$double.eps^0.25)
abs(oneDres$objective-getGlobalOpt("AluffiPentini"))

# Four dimensional optimization task
getDefaultBounds("CosMix4")
oneDres <- optimize1dMulti (f_input=function (x) 
goTest (x, fnName="CosMix4", checkDim=TRUE), 
maxRuns=3, repetitions=5, interval_matrix=matrix(rep(c(-2, 1), 4), ncol=4), 
tol_input=.Machine$double.eps^0.25)
abs(oneDres$objective-getGlobalOpt("CosMix4"))

# Ten dimensional optimization task
getDefaultBounds("Rastrigin")
oneDres <- optimize1dMulti (f_input=function (x) 
goTest (x, fnName="Rastrigin", checkDim=TRUE), 
maxRuns=3, repetitions=5, interval_matrix=matrix(rep(c(-525, 512), 10), ncol=10),
tol_input=.Machine$double.eps^0.25)
abs(oneDres$objective-getGlobalOpt("Rastrigin"))

# Ten dimensional optimization task with higher accuracy
getDefaultBounds("Rastrigin")
oneDres <- optimize1dMulti (f_input=function (x) 
goTest (x, fnName="Rastrigin", checkDim=TRUE), 
maxRuns=3, repetitions=5, interval_matrix=matrix(rep(c(-525, 512), 10), ncol=10),
tol_input=.Machine$double.eps^0.5)
abs(oneDres$objective-getGlobalOpt("Rastrigin"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line