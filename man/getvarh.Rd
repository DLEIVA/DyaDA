\name{getvarh}

\alias{getvarh}

\title{Variance for Landau's Linearity Statistic -h-}

\description{
Function to obtain variance for Landau's h statistic.}

\usage{
   getvarh(X)
}

\arguments{
  \item{X}{Original sociomatrix.}

}

\details{
\code{getvarh} computes variance for Landau's h with the formula \eqn{18 \over{N^3}}.
}

\value{

\item{varianceh}{Variance for Landau's h statistic.}
}

\references{

Landau, H. G. (1951). On dominance relations and the structure of animal societies: I. Effect of inherent characteristics. \emph{Bulletin of Mathematical Biophysics}, \emph{13}, 1-19.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.
}

\seealso{ \code{\link{getexpech}}, \code{\link{getimplandau}}, \code{\link{getlandau}}, \code{\link{linear.hierarchy.test}}.
}

\examples{

  set.seed(123)
  X <- matrix(floor(runif(20*20,min=0,max=20)),nrow=20,byrow=TRUE)
  diag(X)<-0
  getvarh(X)
}

\keyword{misc}

\keyword{htest}

