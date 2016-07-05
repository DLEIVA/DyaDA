\name{getdc}

\alias{getdc}

\title{Directional consistency index - DC -}

\description{
A function to compute the directional consistency index.}


\usage{
   getdc(X)
}

\arguments{
  \item{X}{Original sociomatrix.}
}

\details{
\code{getdc} is only applied for square matrices in which the set of \emph{n} actors is also the set of \emph{n} partners. 
}

\value{
\item{dc}{ In the computation of the directional consistency index the dyadic discrepancies are taken into account. It ranges from 0 to 1 and it can be used as a social reciprocity measure.}
}

\references{

van Hooff, J. A. R. A. M., & Wensing, J. A. B. (1987). Dominance and its behavioral measures in a captive wolf pack. In H. W. Frank (Ed.), \emph{Man and wolf} (pp. 219-252). Dordrecht: Junk Publishers.
}

\author{
David Leiva \email{dleivaur@ub.edu}.
}

\seealso{ \code{\link{getphi}},\code{\link{reciptest1}}
} 

\examples{
##
  X=matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,ncol=3)
}

\keyword{misc}

\keyword{htest}

