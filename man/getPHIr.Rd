\name{getPHIr}
\alias{getPHIr}
\title{Overall asymmetry index -PHIr-}
\description{
Function to obtain the PHIr index by means of the observed sociomatrix.}
\usage{
   getPHIr(X)
}
\arguments{
  \item{X}{Original sociomatrix.}
}
\details{
\code{getPHIr} computes the overall asymmetry index as the sum of the dyadic directional contributions to asymmetry in the whole group, \emph{phirmat}. This measure is a social reciprocity index given that takes into account what is given and what is received by all individuals in the group.
}
\value{
\item{PHIr}{Overall asymmetry index.}
}
\references{
Solanas, A., Leiva, D., Sierra, V., & Salafranca, Ll. (2009). Measuring and making decisions for social reciprocity. \emph{Behavior Research Methods}. \emph{41}, 742-754.
}
\author{
David Leiva \email{dleivaur@ub.edu} & Antonio Solanas \email{antonio.solanas@ub.edu}.
}
\seealso{ \code{\link{getphirmat}}, \code{\link{getPHIrexpec}}, \code{\link{getPHIrSE}}.
} 
\examples{
  X=matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,ncol=3)
}
\keyword{misc}
\keyword{htest}
