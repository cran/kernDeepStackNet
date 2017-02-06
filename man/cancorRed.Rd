\name{cancorRed}
\alias{cancorRed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate first canonical correlation
}
\description{
Calculates the first canonical correlation between two sets of variables. 
}
\usage{
cancorRed(x, y, xcenter = TRUE, ycenter = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
Numeric matrix (n * p1), containing the x coordinates.
}
\item{y}{
Numeric matrix (n * p2), containing the y coordinates.
}
\item{xcenter}{
Logical or numeric vector of length p1, describing any centering to be done on the x values before the analysis. If the argument is TRUE (default), then all columns are centered. Otherwise not adjustment is made.
}
\item{ycenter}{
analogous to xcenter, but for the y values.
}
}
\details{
The canonical correlation analysis seeks linear combinations of the y variables which are well explained by linear combinations of the x variables. The relationship is symmetric as 'well explained' is measured by correlations.
}

\value{
First canonical correlation (numeric vector).
}

\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}
}

\note{
This function is a reduced version of the original Code in the base package \code{\link[stats]{cancor}}. In comparison to the original, it does only calculate the first canonical correlation, but the runtime is reduced.
}

\seealso{
\code{\link[stats]{cancor}}, \code{\link{rdcPart}}, \code{\link{rdcVarOrder}}
}

%\examples{
%
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
