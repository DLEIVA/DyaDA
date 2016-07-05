\name{getepsilon}

\alias{getepsilon}

\title{Generalized reciprocity index - epsilon -}

\description{
Function to obtain the generalized reciprocity index by means of the observed sociomatrix.
}


\usage{
   getepsilon(X)
}

\arguments{
  \item{X}{Original sociomatrix}
}

\details{
\code{getepsilon} is only applied for square matrices in which the set of \emph{n} actors is also the set of \emph{n} partners. 
}

\value{
\item{epsilon}{ The generalized reciprocity index is computed by means of the dyadic and individual contributions to the overall skew-symmetry. It ranges from 0 to 1.}
}

\references{

Solanas, A., Salafranca, Ll., Riba, C., Sierra, V., & Leiva, D. (2006). Quantifying social asymmetric structures. \emph{Behavior Research Methods}, \emph{38}, 390-399. 
}

\author{
David Leiva \email{dleivaur@ub.edu}.
}

\seealso{ \code{\link{getnuj}}, \code{\link{getrationu}}, \code{\link{reciptest1}}
} 

\examples{
##
  X=matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,ncol=3)
}

\keyword{misc}

\keyword{htest}
