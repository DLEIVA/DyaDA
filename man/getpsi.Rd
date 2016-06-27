\name{getpsi}

\alias{getpsi}

\title{Symmetry index - psi -}

\description{
Function to obtain the symmetry index by means of the observed sociomatrix.
}


\usage{
   getpsi(X)
}

\arguments{
  \item{X}{Original sociomatrix}
}

\details{
\code{getpsi} is only applied for square matrices in which the set of \emph{n} actors is also the set of \emph{n} partners. 
}

\value{
\item{psi}{ The symmetry index is based on dyadic discrepancies. The index is the ratio between the sum of squared values due to the symmetrical part and the total sum of squared values.It ranges from 0.5 to 1 and it is a measure of social reciprocity.}
}

\references{

Solanas, A., Salafranca, Ll., Riba, C., Sierra, V., & Leiva, D. (2006). Quantifying social asymmetric structures. \emph{Behavior Research Methods}, \emph{38}, 390-399. 
}

\author{
David Leiva \email{dleivaur@ub.edu}.
}

\seealso{ \code{\link{getphi}},\code{\link{reciptest}}
} 

\examples{
##
  X=matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,ncol=3)
}

\keyword{misc}

\keyword{htest}
