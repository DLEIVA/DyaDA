\name{getphirmat}
\alias{getphirmat}
\title{Dyadic directional asymmetry matrix -phirmat-}
\description{
Function to obtain the matrix of dyadic directional contributions to asymmetry measures by means of the observed sociomatrix.}
\usage{
   getphirmat(X)
}
\arguments{
  \item{X}{Original sociomatrix.}
}
\details{
\code{getphirmat} computes the dyadic directional contributions to asymmetry. The \emph{phirmat} matrix is required in order to obtain the other group, dyadic and individual measures of social reciprocity here proposed.
}
\value{
\item{phirmat}{Dyadic directional asymmetry matrix.}
}
\references{
Solanas, A., Leiva, D., Sierra, V., & Salafranca, Ll. (2009). Measuring and making decisions for social reciprocity. \emph{Behavior Research Methods}. \emph{41}, 742-754.
}
\author{
David Leiva \email{dleivaur@ub.edu} & Antonio Solanas \email{antonio.solanas@ub.edu}.
}
\seealso{ \code{\link{getPHIr}}, \code{\link{getphii}}, \code{\link{getphij}}, \code{\link{getPHIij}}.
} 
\examples{
  X=matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,ncol=3)
}
\keyword{misc}
\keyword{htest}
