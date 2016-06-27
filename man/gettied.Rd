\name{gettied}

\alias{gettied}

\title{Dyads with tied relationships}

\description{
Function to obtain the number of tied relationships of the observed sociomatrix.}

\usage{
   gettied(X,relfreq=TRUE)
}

\arguments{
  \item{X}{Original sociomatrix.}
  \item{relfreq}{Logical. If TRUE the function computes percentage of tied relationships.}

}

\details{
\code{gettied} computes the number of dyads in which the individuals receive and address an equal number of interactions.
}

\value{

\item{tied.dyads}{Number of dyads with tied relationships.}
\item{perc.tied.dyads}{Percentage of dyads with tied relationships.}
}

\references{

van Hooff, J. A. R. A. M., & Wensing, J. A. B. (1987). Dominance and its behavioral measures in a captive wolf pack. In H. W. Frank (Ed.), \emph{Man and wolf} (pp. 219-252). Dordrecht: Junk Publishers.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.
}

\seealso{ \code{\link{getonewayrel}}, \code{\link{gettwowayrel}}, \code{\link{getunknown}}.
}

\examples{

  set.seed(123)
  X <- matrix(floor(runif(20*20,min=0,max=20)),nrow=20,byrow=TRUE)
  diag(X)<-0
  gettied(X)
}

\keyword{misc}

\keyword{htest}

