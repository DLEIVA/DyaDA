\name{getlandau}

\alias{getlandau}

\title{Landau's Linearity Index -h-}

\description{
Function to obtain Landau's h index for quantifying linearity in social dominance hierarchies.}

\usage{
   getlandau(X)
}

\arguments{
  \item{X}{Original sociomatrix.}

}

\details{
\code{getlandau} computes Landau's h index (1951) for quantifying the degree of linearity in social dominance hierarchies. The linearity index is defined as:
 \tabular{ll}{
	\tab \eqn{
  		h = {12  \over {N^3-N}} \sum\limits_{i = 1}^N {\left(V_i-{{N-1} \over 2} \right)^2}}, \cr}

where \eqn{V_i} is the number of individuals dominated by the individual \emph{i}. The elements \eqn{V_i} of the vector \emph{V} is just the row sum of matrix of abilities, which is computed by the function called \emph{freq2abil}.

}

\value{

\item{landauh}{Landau's index h.}
}

\references{

Landau, H. G. (1951). On dominance relations and the structure of animal societies: I. Effect of inherent characteristics. \emph{Bulletin of Mathematical Biophysics}, \emph{13}, 1-19.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.
}

\seealso{\code{\link{freq2abil}}, \code{\link{getexpech}}, \code{\link{getimplandau}}, \code{\link{getvarh}}, \code{\link{linear.hierarchy.test}}.
}

\examples{

  set.seed(123)
  X <- matrix(floor(runif(20*20,min=0,max=20)),nrow=20,byrow=TRUE)
  diag(X)<-0
  getlandau(X)
}

\keyword{misc}

\keyword{htest}

