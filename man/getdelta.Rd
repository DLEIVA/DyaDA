\name{getdelta}

\alias{getdelta}

\title{Delta index - delta -}

\description{
Function to obtain the delta index by means of the observed sociomatrix.
}


\usage{
   getdelta(X)
}

\arguments{
  \item{X}{Original sociomatrix}
}

\details{
\code{getdelta} is only applied for square matrices in which the set of \emph{n} actors is also the set of \emph{n} partners. 
}

\value{
\item{delta}{ The delta index is the ratio of overall skew-symmetry and symmetry indices, phi and psi, respectively. It ranges from 0 to 1.}
}

\references{

Solanas, A., Salafranca, Ll., Riba, C., Sierra, V., & Leiva, D. (2006). Quantifying social asymmetric structures. \emph{Behavior Research Methods}, \emph{38}, 390-399. 
}

\author{
David Leiva \email{dleivaur@ub.edu}.
}

\seealso{ \code{\link{getphi}}, \code{\link{getpsi}}, \code{\link{reciptest}}
} 

\examples{
##
  X=matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,ncol=3)
}

\keyword{misc}

\keyword{htest}
