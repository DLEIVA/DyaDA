\name{getPHIrexpec}
\alias{getPHIrexpec}
\title{Mathematical expectancy of the PHIr statistic -PHIrexpec-}
\description{
Function to obtain mathematical expectancy of the PHIr statistic under a specified null hypothesis regarding social reciprocity and by means of the observed sociomatrix.}
\usage{
   getPHIrexpec(X,pi)
}
\arguments{
  \item{X}{Original sociomatrix.}
  \item{pi}{Matrix of probabilities \emph{Pij} that specifies the null hypothesis regarding social reciprocity.}  
}
\details{
\code{getPHIrexpec} computes the mathematical expectancy of the PHIr statistic under any null hypothesis. The null hypothesis of complete reciprocation could be of interest to social researchers. This hypothesis can be represented as \eqn{Pij = .5 for all i and j}.
}
\value{
\item{PHIrexpec}{Mathematical expectancy of the PHIr statistic under the specified null hypothesis.}
}
\references{
Solanas, A., Leiva, D., Sierra, V., & Salafranca, Ll. (2009). Measuring and making decisions for social reciprocity. \emph{Behavior Research Methods}. \emph{41}, 742-754.
}
\author{
David Leiva \email{dleivaur@ub.edu} & Antonio Solanas \email{antonio.solanas@ub.edu}.
}
\seealso{ \code{\link{getPHIr}}.
} 
\examples{
  X=matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,ncol=3)
  pi=matrix(c(0, 0.5, 0.5, 0.5, 0, 0.5, 0.5, 0.5, 0),nrow=3,ncol=3)
}
\keyword{misc}
\keyword{htest}
