\name{crossprodRcpp}
\alias{crossprodRcpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculates the cross product of a matrix
}
\description{
This is an computational efficient help function to calculate the GCV loss of the KDSN. Calculates the matrix product A^T A. It is written in C and is usually faster than the regular R version.
}
\usage{
  crossprodRcpp(A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
    Numeric Matrix. Integer matrices are not allowed.
  }
}
%\details{
%}
\value{
List with the numeric matrix A^T A as element.
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
\code{\link[base]{crossprod}}
}

\examples{
A <- matrix(seq(-1, 1, length.out=9), nrow=3, ncol=3)
crossprodRcpp(A)[[1]]
all.equal(crossprodRcpp(A)[[1]], crossprod(A))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models & regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line