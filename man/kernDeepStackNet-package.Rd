\name{kernDeepStackNet-package}
\alias{kernDeepStackNet-package}
\alias{kernDeepStackNet}
\docType{package}
\title{
Kernel deep stacking networks with random Fourier transformation
}
\description{
Contains functions for estimation and model selection of kernel deep stacking networks (KDSNs). KDSNs are a supervised learning method, which can be used for continuous or binary responses. The model specification follows the approach of Huang et al. (2013), which is based on a series of kernel ridge regression models to random Fourier transformed input data. 

The model selection includes model-based optimization of arbitrary loss functions. All help functions are also available for customized modeling, but it is recommended to use the higher level functions. The main functions are 
\itemize{
%\item{\code{\link{tuneMboKDSN}}: }{Model selection of tuning parameters of 
%kernel deep stacking networks with model-based optimization of generalized cross-valdiation loss function.}
%\item{\code{\link{tuneMboCvKDSN}}: }{Model selection of tuning parameters of 
%kernel deep stacking networks with model-based optimization using cross-valdiation with arbitrary loss functions.}
\item{\code{\link{fitKDSN}}: }{Fits kernel deep stacking networks.}
\item{\code{\link{tuneMboLevelCvKDSN}}: }{Selection of tuning parameters of 
kernel deep stacking networks with model-based optimization using cross-valdiation, arbitrary loss functions and pre-specified number of levels. All other tuning parameters are free to vary across the levels.}
\item{\code{\link{tuneMboLevelGcvKDSN}}: }{Selection of tuning parameters of 
kernel deep stacking networks with model-based optimization using generalized cross validation score, arbitrary loss functions and pre-specified number of levels. All other tuning parameters are free to vary across the levels.}
\item{\code{\link{tuneMboSharedCvKDSN}}: }{Selection of tuning parameters of 
kernel deep stacking networks with model-based optimization using cross validation or test set and arbitrary loss functions. The number of levels is included in tuning. All other tuning parameters are shared across the levels.}
\item{\code{\link{tuneMboSharedSubsetKDSN}}: }{Selection of tuning parameters of 
kernel deep stacking networks with model-based optimization using test sets and arbitrary loss functions. The number of levels is included in tuning. All other tuning parameters are shared across the levels. In contrast to function \code{\link{tuneMboSharedCvKDSN}} an ensemble of (sparse) KDSNs is estimated.}
%\item{\code{\link{fineTuneCvKDSN}}: }{Finds better suited random Fourier transformation weights, regarding arbitrary cross-validated loss functions, by random search given an estimated KDSN structure.}
}
For examples and further information, see the corresponding readme pages of the main functions.
}
\details{
\tabular{ll}{
Package: \tab kernDeepStackNet \cr
Type: \tab Package \cr
Version: \tab 2.0.2 \cr
Date: \tab 2017-05-31 \cr
License: \tab GPL-3 \cr
}
%%~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Maintainer: Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}

Matthias Schmid \email{matthias.schmid@imbie.uni-bonn.de}
}
\references{

Thomas Welchowski, Matthias Schmid, (2016),
\emph{A framework for parameter estimation and model selection in kernel deep stacking networks},
Artificial Intelligence in Medicine, volume 70, pages 31-40

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

Victor Picheny, David Ginsbourger, Yann Richet, (2012),
\emph{Quantile-based optimization of Noisy Computer Experiments with Tunable Precision},
HAL-archives-ouvertes.fr, hal-00578550v3

Nitish Srivastava, Geoffrey Hinton, Alex Krizhevsky et al., (2014),
\emph{Dropout: A simple way to prevent neural networks from overfitting},
Journal of Machine Learning Research, volume 15, pages 1929-1958

David Lopez-Paz and Philipp Hennig and Bernhard Sch\"olkopf, (2013),
\emph{The randomized dependence coefficient},
https://arxiv.org/, reference arXiv:1304.7717

Luca Scrucca, (2013), 
\emph{GA: A Package for Genetic Algorithms in R},
Journal of Statistical Software, issue 4, volume 53, pages 1-37

Constantino Tsallis, Daniel A. Stariolo, (1996),
\emph{Generalized simulated annealing},
Elsevier, Physica A, volume 233, pages 395-406

}

\keyword{package}

%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R ~~
%%~~ documentation directory ~~
%\seealso{
%%~~ Optional links to other man pages, e.g. ~~
%%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{
%
%}
