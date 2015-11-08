\name{getEigenValuesRcpp}
\alias{getEigenValuesRcpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the eigenvalues of a matrix
}
\description{
This is an computational efficient help function to calculate the GCV loss of the KDSN. It is written in C and is faster than the regular R version.
}
\usage{
getEigenValuesRcpp(M)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{M}{
Numeric Matrix.
}
}
%\details{
%}
\value{
Numeric vector with increasing eigenvalues.
}
%\references{
%}
\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}
}

%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{eigen}}
}

%\examples{
%
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line