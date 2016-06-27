\name{getimplandau}

\alias{getimplandau}

\title{Improved Landau's Linearity Index -h'-}

\description{
Function to obtain improved Landau's h' index for quantifying the degree of linearity in social dominance hierarchies proposed by de Vries (1995).}

\usage{
   getimplandau(X)
}

\arguments{
  \item{X}{Original sociomatrix.}

}

\details{
\code{getimplandau} computes improved Landau's h' index in dominance hierarchies containing unknown or tied relationships. The index proposed by de Vries (1995) can be computed as follows:

 \tabular{ll}{
	\tab \eqn{
  		h' = h + {6  \over {N^3-N}} \times u}, \cr
}

where \emph{h} is Landau's index and \emph{u} is the number of unknown dyads. 

}

\value{

\item{improvedlandauh}{Improved Landau's index h'.}
}

\references{

de Vries, H. (1995). An improved test of linearity in dominance hierarchies containing unknown or tied relationships. \emph{Animal Behaviour}, \emph{50}, 1375-1389.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.
}

\seealso{ \code{\link{getexpech}}, \code{\link{getlandau}}, \code{\link{getvarh}}, \code{\link{linear.hierarchy.test}}.
}

\examples{

  set.seed(123)
  X <- matrix(floor(runif(20*20,min=0,max=20)),nrow=20,byrow=TRUE)
  diag(X)<-0
  getimplandau(X)
}

\keyword{misc}

\keyword{htest}

