\name{getnu}

\alias{getnu}

\title{Matrix of unweighted contributions to skew-symmetry - nu -}

\description{
Function to obtain the matrix of unweighted contributions to asymmetry by means of the observed sociomatrix.
}


\usage{
   getnu(X)
}

\arguments{
  \item{X}{Original sociomatrix}
}

\details{
\code{getnu} is only applied for square matrices in which the set of \emph{n} actors is also the set of \emph{n} partners. 
}

\value{
\item{nu}{ An \emph{n} x \emph{n} matrix of unweighted contributions to asymmetry is provided. By means of this matrix dyadic and individual contributions to skew-symmetry weighted by interaction frequency will be computed.}
}

\references{

Solanas, A., Salafranca, Ll., Riba, C., Sierra, V., & Leiva, D. (2006). Quantifying social asymmetric structures. \emph{Behavior Research Methods}, \emph{38}, 390-399. 
}

\author{
David Leiva \email{dleivaur@ub.edu}.
}

\seealso{ \code{\link{getnuj}},\code{\link{getrationu}},\code{\link{reciptest1}}
} 

\examples{
##
  X=matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,ncol=3)
}

\keyword{misc}

\keyword{htest}
