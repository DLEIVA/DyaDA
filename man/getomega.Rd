\name{getomega}

\alias{getomega}

\title{Matrix of dyadic balanced reciprocity - omega -}

\description{
Function to obtain the matrix of dyadic balanced reciprocity by means of the observed sociomatrix.
}


\usage{
   getomega(X)
}

\arguments{
  \item{X}{Original sociomatrix}
}

\details{
\code{getomega} is only applied for square matrices in which the set of \emph{n} actors is also the set of \emph{n} partners. 
}

\value{
\item{omega}{ An \emph{n} x \emph{n} matrix of dyadic balanced reciprocity is provided. The elements of the matrix are obtained by dividing dyadic contributions to skew-symmetry by dyadic contributions to symmetry, rationu and ratiolambda, respectively.}
}

\references{

Solanas, A., Salafranca, Ll., Riba, C., Sierra, V., & Leiva, D. (2006). Quantifying social asymmetric structures. \emph{Behavior Research Methods}, \emph{38}, 390-399. 
}

\author{
David Leiva \email{dleivaur@ub.edu}.
}

\seealso{ \code{\link{getratiolambda}},\code{\link{getrationu}},\code{\link{reciptest1}}
} 

\examples{
##
  X=matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,ncol=3)
}

\keyword{misc}

\keyword{htest}
