\name{getPHIij}
\alias{getPHIij}
\title{Dyadic contributions to asymmetry matrix -PHIij-}
\description{
Function to obtain the matrix of dyadic contributions to asymmetry measures by means of the observed sociomatrix.
}
\usage{
   getPHIij(X)
}
\arguments{
  \item{X}{Original sociomatrix.}
}
\details{
\code{getPHIij} computes the dyadic contributions to asymmetry. The \emph{phirmat} matrix is required to obtain this symmetrical matrix.
}
\value{
\item{PHIij}{Dyadic contributions to asymmetry matrix.}
}
\references{
Solanas, A., Leiva, D., Sierra, V., & Salafranca, Ll. (2009). Measuring and making decisions for social reciprocity. \emph{Behavior Research Methods}. \emph{41}, 742-754.
}
\author{
David Leiva \email{dleivaur@ub.edu} & Antonio Solanas \email{antonio.solanas@ub.edu}.
}
\seealso{ \code{\link{getphirmat}}.
} 
\examples{
  X=matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,ncol=3)
}
\keyword{misc}
\keyword{htest}
