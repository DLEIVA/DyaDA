\name{getexpech}

\alias{getexpech}

\title{Mathematical Expectancy for Landau's Linearity Index -h-}

\description{
Function to obtain mathematical expectancy for Landau's h index.}

\usage{
   getexpech(X)
}

\arguments{
  \item{X}{Original sociomatrix.}

}

\details{
\code{getexpech} computes mathematical expectancy for Landau's h with the formula \eqn{3 \over {N+1}}.
}

\value{

\item{expectedh}{Mathematical expectancy for Landau's index h.}
}

\references{

Landau, H. G. (1951). On dominance relations and the structure of animal societies: I. Effect of inherent characteristics. \emph{Bulletin of Mathematical Biophysics}, \emph{13}, 1-19.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.
}

\seealso{ \code{\link{getimplandau}}, \code{\link{getlandau}}, \code{\link{getvarh}}, \code{\link{linear.hierarchy.test}}.
}

\examples{

  set.seed(123)
  X <- matrix(floor(runif(20*20,min=0,max=20)),nrow=20,byrow=TRUE)
  diag(X)<-0
  getexpech(X)
}

\keyword{misc}

\keyword{htest}

