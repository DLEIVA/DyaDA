\name{getonewayrel}

\alias{getonewayrel}

\title{Dyads with one way relationships}

\description{
Function to obtain the number of dyads with one way relationships of the observed sociomatrix.}

\usage{
   getonewayrel(X, relfreq=TRUE)
}

\arguments{
  \item{X}{Original sociomatrix.}
  \item{relfreq}{Logical. If TRUE the function computes percentage of dyads with one way relationships.}
}

\details{
\code{getonewayrel} computes the number of dyads in which one member addresses behaviours and the other does not.
}

\value{

\item{oneway.dyads}{Number of dyads with one way relationships.}
\item{perc.oneway.dyads}{Percentage of dyads with one way relationships.}
}

\references{

van Hooff, J. A. R. A. M., & Wensing, J. A. B. (1987). Dominance and its behavioral measures in a captive wolf pack. In H. W. Frank (Ed.), \emph{Man and wolf} (pp. 219-252). Dordrecht: Junk Publishers.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.
}

\seealso{ \code{\link{gettied}}, \code{\link{gettwowayrel}}, \code{\link{getunknown}}.
}

\examples{

  set.seed(123)
  X <- matrix(floor(runif(20*20,min=0,max=20)),nrow=20,byrow=TRUE)
  diag(X)<-0
  getonewayrel(X)
}

\keyword{misc}

\keyword{htest}

