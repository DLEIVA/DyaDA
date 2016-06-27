\name{gettwowayrel}

\alias{gettwowayrel}

\title{Dyads with two way relationships}

\description{
Function to obtain the number of dyads with two way relationships of the observed sociomatrix.}

\usage{
   gettwowayrel(X, relfreq=TRUE)
}

\arguments{
  \item{X}{Original sociomatrix.}
  \item{relfreq}{Logical. If TRUE the function computes relative frequencies of two way relationships.}
}

\details{
\code{gettwowayrel} computes the number of dyads in which both members address and receive behaviours.
}

\value{

\item{twoway.dyads}{Number of dyads with two way relationships.}
\item{perc.twoway.dyads}{Percentage of dyads with two relationships.}
}

\references{

van Hooff, J. A. R. A. M., & Wensing, J. A. B. (1987). Dominance and its behavioral measures in a captive wolf pack. In H. W. Frank (Ed.), \emph{Man and wolf} (pp. 219-252). Dordrecht: Junk Publishers.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.
}

\seealso{ \code{\link{getonewayrel}}, \code{\link{gettied}}, \code{\link{getunknown}}.
}

\examples{

  set.seed(123)
  X <- matrix(floor(runif(20*20,min=0,max=20)),nrow=20,byrow=TRUE)
  diag(X)<-0
  gettwowayrel(X)
}

\keyword{misc}

\keyword{htest}

