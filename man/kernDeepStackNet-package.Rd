\name{kernDeepStackNet-package}
\alias{kernDeepStackNet-package}
\alias{kernDeepStackNet}
\docType{package}
\title{
Kernel deep stacking networks with random Fourier transformation
}
\description{
Contains functions for estimation and model selection of kernel deep stacking networks for supervised learning. The model specification follows the approach of Huang et al. (2013), which is based on a series of kernel ridge regression models to random Fourier transformed input data.

The model selection includes direct optimization or model-based alternatives of arbitrary loss functions. All help functions are also available for customized modeling, but it is recommended to use the higher level functions. The main functions are 
\itemize{
\item{\code{\link{fitKDSN}}: }{Fits kernel deep stacking networks.}
\item{\code{\link{tuneMboKDSN}}: }{Model selection of tuning parameters of 
kernel deep stacking networks with model-based optimization of generalized cross-valdiation loss function.}
\item{\code{\link{tuneMboCvKDSN}}: }{Model selection of tuning parameters of 
kernel deep stacking networks with model-based optimization using cross-valdiation with arbitrary loss functions.}
\item{\code{\link{tuneMboLevelCvKDSN}}: }{Model selection of tuning parameters of 
kernel deep stacking networks with model-based optimization using cross-valdiation, arbitrary loss functions and pre-specified number of levels.}
\item{\code{\link{fineTuneCvKDSN}}: }{Finds better suited random Fourier transformation weights, regarding arbitrary cross-validated loss functions, by random search given an estimated KDSN structure.}
}
}
\details{
\tabular{ll}{
Package: \tab kernDeepStackNet\cr
Type: \tab Package\cr
Version: \tab 1.0.1\cr
Date: \tab 2016-02-01 \cr
License: \tab GPL-3 \cr
}
%%~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Maintainer: Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}

Matthias Schmid \email{matthias.schmid@imbie.uni-bonn.de}
}
\references{
Po-Seng Huang and Li Deng and Mark Hasegawa-Johnson and Xiaodong He, (2013), 
\emph{Random Features for kernel deep convex network},
Proceedings IEEE International Conference on Acoustics, Speech, and 
Signal Processing (ICASSP)

Simon N. Wood, (2006), 
\emph{Generalized Additive Models: An Introduction with R},
Taylor & Francis Group LLC

R. Brent, (1973), 
\emph{Algorithms for Minimization without Derivatives},
Englewood Cliffs N.J.: Prentice-Hall

Donald R. Jones and Matthias Schonlau and William J. Welch, (1998), 
\emph{Efficient Global Optimization of Expensive Black-Box Functions},
Journal of Global Optimization 13: pages 455-492

Krige DG, (1951),
\emph{A Statistical Approach to Some Basic Mine Valuation Problems 
on the Witwatersrand},
Journal of the Chemical, Metallurgical and Mining Society of South Africa, 
52(6), 119-139

Olivier Roustant and David Ginsbourger and Yves Deville, (2012), 
\emph{DiceKriging, DiceOptim: Two R Packages for the Analysis of Computer Experiments 
by Kriging-Based Metamodeling and Optimization},
Journal of Statistical Software, Volume 51, Issue 1

Michael Stein, (1987),
\emph{Large Sample Properties of Simulations Using Latin Hypercube Sampling},
Technometrics. 29, 143-151

Carl Edward Rasmussen and Christopher K. I. Williams, (2006), 
\emph{Gaussian Processes for Machine Learning},
Massachusetts Institute of Technology

Jerome Friedman and Trevor Hastie and Rob Tibshirani, (2008), 
\emph{Regularization Paths for Generalized Linear Models via Coordinate Descent},
Department of Statistics, Stanford University
}
\keyword{ package }

%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R ~~
%%~~ documentation directory ~~
%\seealso{
%%~~ Optional links to other man pages, e.g. ~~
%%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%%\examples{
%%~~ simple examples of the most important functions ~~
%%}