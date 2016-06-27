\name{getcirc}

\alias{getcirc}

\title{Number of circular triads}

\description{
Function to obtain the number of circular triads of an observed matrix.}

\usage{
   getcirc(X)
}

\arguments{
  \item{X}{Original sociomatrix.}

}

\details{
\code{getcirc} computes the number of circular triads in a sociomatrix with the following formula:
 \tabular{ll}{
	\tab \eqn{
  		d = {{N \left( N-1 \right) \left( 2N-1 \right)}  \over 12} - {{1\over2} \sum\limits_{i = 1}^N {S_i}^2}}, \cr
} 
where \eqn{S_i} is the row sum of the matrix of dominances.

This function it is also used in \emph{getkendall} in order to compute the circularity index proposed by Kendall and Babington Smith (1940).
}

\value{

\item{circular}{Number of circular triads.}
}

\references{

Kendall, M. G., & Babington Smith, B. (1940). On the method of paired comparisons. \emph{Biometrika}, \emph{31}, 324-345.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.
}

\seealso{ \code{\link{getexpeccirc}}, \code{\link{getkendall}}, \code{\link{getpvkendall}}.
}

\examples{

  set.seed(123)
  X <- matrix(floor(runif(20*20,min=0,max=20)),nrow=20,byrow=TRUE)
  diag(X)<-0
  getcirc(X)
}

\keyword{misc}

\keyword{htest}

