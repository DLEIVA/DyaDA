\name{getmaxcirc}

\alias{getmaxcirc}

\title{Maximum number of circular triads}

\description{
Function to obtain the maximum number of circular triads of an observed matrix.}

\usage{
   getmaxcirc(X)
}

\arguments{
  \item{X}{Original sociomatrix.}

}

\details{
\code{getmaxcirc} computes the maximum number of circular triads in a sociomatrix following original formulae proposed by Kendall and Babington Smith (1940):
 \tabular{ll}{
	\tab \eqn{
  		maxd = {{1 \over 24} \left( N^3 - N \right)}} {for N odd}, \cr
        \tab \eqn{
                maxd = {{1 \over 24} \left( N^3 - 4N \right)}} {for N even}. \cr
}

This function is used in \emph{getkendall} in order to compute the K index for quantifying the degree of linearity in social dominance hierarchies.
}

\value{

\item{maxcirc}{Maximum number of circular triads.}
}

\references{

Kendall, M. G., & Babington Smith, B. (1940). On the method of paired comparisons. \emph{Biometrika}, \emph{31}, 324-345.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.
}

\seealso{ \code{\link{getcirc}}, \code{\link{getexpeccirc}}, \code{\link{getkendall}}, \code{\link{getpvkendall}}.
}

\examples{

  set.seed(123)
  X <- matrix(floor(runif(20*20,min=0,max=20)),nrow=20,byrow=TRUE)
  diag(X)<-0
  getmaxcirc(X)
}

\keyword{misc}

\keyword{htest}

