\name{kernDeepStackNet_getEigenValuesRcpp}
\alias{kernDeepStackNet_getEigenValuesRcpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the eigenvalues of a matrix
}
\description{
This is an computational efficient help function to calculate the GCV loss of the KDSN. It is written in C and is faster than the regular R version.
}
\usage{
kernDeepStackNet_getEigenValuesRcpp(M)
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
Numeric vector eigenvalues in increasing order.
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

\examples{
# Define simple binary matrix
M <- cbind(c(1,-1), c(-1,1))

# Calculate eigen values with standard function
eigenVal <- eigen(M, only.values=TRUE)$values

# Calculate eigen values with Rcpp version
eigenValRcpp <- rev(getEigenValuesRcpp(M))

# Check equality
all.equal(eigenVal, eigenValRcpp) # TRUE

# Analyse run time of both variants
if(tryCatch(require("microbenchmark", quietly=TRUE, warn.conflicts=FALSE, 
character.only = TRUE), warning=function(w) FALSE, error=function(e) FALSE)) {
  microbenchmark(Standard=eigen(M, only.values=TRUE), Rcpp=getEigenValuesRcpp(M))
}
# -> Rcpp function is on average about 30 times faster than the standard R function
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line