\name{getunknown}

\alias{getunknown}

\title{Dyads with unknown relationships}

\description{
Function to obtain the number of dyads with unknown relationships of the observed sociomatrix.}

\usage{
   getunknown(X,relfreq=TRUE)
}

\arguments{
  \item{X}{Original sociomatrix.}
  \item{relfreq}{Logical. If TRUE, the function computes percentage of unknown relationships.}
}

\details{
\code{getunknown} computes the number of dyads with 0 interactions.
}

\value{

\item{unknown.dyads}{Number of dyads with unknown relationships.}
\item{perc.unknown.dyads}{Percentage of dyads with unknown relationships.}

}

\references{

van Hooff, J. A. R. A. M., & Wensing, J. A. B. (1987). Dominance and its behavioral measures in a captive wolf pack. In H. W. Frank (Ed.), \emph{Man and wolf} (pp. 219-252). Dordrecht: Junk Publishers.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.
}

\seealso{ \code{\link{getonewayrel}}, \code{\link{gettied}}, \code{\link{gettwowayrel}}.
}

\examples{

  set.seed(123)
  X <- matrix(floor(runif(20*20,min=0,max=20)),nrow=20,byrow=TRUE)
  diag(X)<-0
  getunknown(X)
}

\keyword{misc}

\keyword{htest}

