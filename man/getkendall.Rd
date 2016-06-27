\name{getkendall}

\alias{getkendall}

\title{Kendall's Linearity Index}

\description{
Function to obtain Kendall's linearity index for measuring linear hierarchies.}

\usage{
   getkendall(X)
}

\arguments{
  \item{X}{Original sociomatrix.}

}

\details{
\code{getkendall} computes Kendall's linearity index as follows:
 \tabular{ll}{
	\tab \eqn{
  		K = {1 - {d \over maxd}}}, \cr}
where d is the number of circular triads in the original sociomatrix and maxd is the maximum number of circular triads.
}

\value{

\item{K}{Kendall's linearity index.}
}

\references{

Kendall, M. G., & Babington Smith, B. (1940). On the method of paired comparisons. \emph{Biometrika}, \emph{31}, 324-345.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.
}

\seealso{ \code{\link{getcirc}}, \code{\link{getexpeccirc}}, \code{\link{getmaxcirc}}, \code{\link{getpvkendall}}.
}

\examples{

  set.seed(123)
  X <- matrix(floor(runif(20*20,min=0,max=20)),nrow=20,byrow=TRUE)
  diag(X)<-0
  getkendall(X)
}

\keyword{misc}

\keyword{htest}

