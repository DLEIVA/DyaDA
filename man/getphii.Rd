\name{getphii}
\alias{getphii}
\title{Individuals' contributions to asymmetry -phii-}
\description{
Function to obtain individuals' contributions as actors by means of the observed sociomatrix.}
\usage{
getphii(X)
}
\arguments{
\item{X}{Original sociomatrix.}
}
\details{
\code{getphii} computes individuals' contributions to asymmetry as actors.
}
\value{
\item{phii}{Individuals' contributions to asymmetry as actors.}
}
\references{
Solanas, A., Leiva, D., Sierra, V., & Salafranca, Ll. (2009). Measuring and making decisions for social reciprocity. \emph{Behavior Research Methods}. \emph{41}, 742-754.
}
\author{
David Leiva \email{dleivaur@ub.edu} & Antonio Solanas \email{antonio.solanas@ub.edu}.
}
\seealso{\code{\link{getphirmat}},\code{\link{getphij}}.
}
\examples{
X=matrix(c(0,1,2,1,0,1,3,1,0),nrow=3,ncol=3)
}
\keyword{misc}
\keyword{htest}
