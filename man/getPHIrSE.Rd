\name{getPHIrSE}
\alias{getPHIrSE}
\title{Standard error of the PHIr statistic -PHIrSE-}
\description{
Function to obtain standard error of the PHIr statistic under a specified null hypothesis regarding social reciprocity and by means of the observed sociomatrix.}
\usage{
   getPHIrSE(X,pi)
}
\arguments{
  \item{X}{Original sociomatrix.}
  \item{pi}{Matrix of probabilities \emph{Pij} that specifies the null hypothesis regarding social reciprocity.}
}
\details{
\code{getPHIrSE} computes the standard error of the PHIr statistic under any null hypothesis. Null hypothesis of complete reciprocation could be of interest to social researchers. This hypothesis can be represented as \eqn{Pij = .5 for all i and j}.
}
\value{
\item{PHIrSE}{Standard error of the PHIr statistic under the specified null hypothesis.}
}
\references{
Solanas, A., Leiva, D., Sierra, V., & Salafranca, Ll. (2009). Measuring and making decisions for social reciprocity. \emph{Behavior Research Methods}. \emph{41}, 742-754.
}
\author{
David Leiva \email{dleivaur@ub.edu} & Antonio Solanas \email{antonio.solanas@ub.edu}.
}
\seealso{ \code{\link{getPHIr}}, \code{\link{getPHIrexpec}}.
}
\examples{
  X=matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,ncol=3)
  pi=matrix(c(0, 0.5, 0.5, 0.5, 0, 0.5, 0.5, 0.5, 0),nrow=3,ncol=3)
}
\keyword{misc}
\keyword{htest}
