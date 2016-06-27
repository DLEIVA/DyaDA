\name{getphij}
\alias{getphij}
\title{Individuals' contributions to asymmetry -phij-}
\description{
Function to obtain individuals' contributions as partners by means of the observed sociomatrix.}
\usage{
   getphij(X)
}
\arguments{
  \item{X}{Original sociomatrix.}
}
\details{
\code{getphij} computes individuals' contribution to asymmetry as partners.
}
\value{
\item{phij}{Individuals' contributions to asymmetry as partners.}
}
\references{
Solanas, A., Leiva, D., Sierra, V., & Salafranca, Ll. (2009). Measuring and making decisions for social reciprocity. \emph{Behavior Research Methods}. \emph{41}, 742-754.
}
\author{
David Leiva \email{dleivaur@ub.edu} & Antonio Solanas \email{antonio.solanas@ub.edu}.
}
\seealso{ \code{\link{getphirmat}}, \code{\link{getphii}}.
} 
\examples{
  X=matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,ncol=3)
}
\keyword{misc}
\keyword{htest}
