\name{EImod}
\alias{EImod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Expected improvement criterion replacement function
}
\description{
This function is a replacement of the function \code{\link[DiceOptim]{EI}} in package DiceOptim. 
% The numerical check contained in EI() caused execution errors and is therefore removed.
}
\usage{
EImod (x, model, plugin = NULL, type = "UK", minimization = TRUE, 
    envir = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
a vector representing the input for which one wishes to calculate EI,
}
\item{model}{
an object of class \code{\link[DiceKriging]{km}}
}
\item{plugin}{
optional scalar: if provided, it replaces the minimum of the current observations,
}
\item{type}{
"SK" or "UK" (by default), depending whether uncertainty related to trend estimation has to be taken into account,
}
\item{minimization}{
logical specifying if EI is used in minimiziation or in maximization,
}
\item{envir}{
an optional environment specifying where to assign intermediate values for future gradient calculations. Default is NULL.
}
}
%\details{
%The line with the numerical check was:
%if (kriging.sd/sqrt(model@covariance@sd2) < 1e-06) {
%        res <- 0
%        xcr <- xcr.prob <- xcr.dens <- NULL
%    }
%The check in the if-clause sometimes does not yield either TRUE or FALSE.
%The if-clause is removed and there is no numerical check.
%}

%return(list(par=c(ego_result), value=attr(ego_result, "funcVal")))

\value{
The expected improvement, defined as

\eqn{EI(x) := E[( min(Y(X)) - Y(x))^{+} | Y(X)=y(X)]},

where X is the current design of experiments and Y is the random process assumed to have generated the objective function y. If a plugin is specified, it replaces

\eqn{min(Y(X))}

in the previous formula.
}

%\references{
%}

\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}
}

\note{
For further examples see DiceOptim::EI
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[DiceOptim]{EI}}
}

%\examples{
%
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line