\name{getexpeccirc}

\alias{getexpeccirc}

\title{Expected number of circular triads}

\description{
Function to obtain the expected number of circular triads of an observed matrix.}

\usage{
   getexpeccirc(X)
}

\arguments{
  \item{X}{Original sociomatrix.}

}

\details{
\code{getexpeccirc} computes the expected number of circular triads in a sociomatrix according the formula \eqn{{\left( N \right) \left( N-1 \right) \left( N-2 \right)} \over 24}.
}

\value{

\item{expectedcirc}{Expected number of circular triads.}
}

\references{

Kendall, M. G., & Babington Smith, B. (1940). On the method of paired comparisons. \emph{Biometrika}, \emph{31}, 324-345.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.
}

\seealso{ \code{\link{getcirc}}, \code{\link{getkendall}}, \code{\link{getmaxcirc}}, \code{\link{getpvkendall}}.
}

\examples{

  set.seed(123)
  X <- matrix(floor(runif(20*20,min=0,max=20)),nrow=20,byrow=TRUE)
  diag(X)<-0
  getexpeccirc(X)
}

\keyword{misc}

\keyword{htest}

