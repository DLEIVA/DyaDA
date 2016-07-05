\name{getlambda}

\alias{getlambda}

\title{Matrix of unweighted contributions to symmetry - lambda -}

\description{
Function to obtain the matrix of unweighted contributions to symmetry by means of the observed sociomatrix.
}


\usage{
   getlambda(X)
}

\arguments{
  \item{X}{Original sociomatrix}
}

\details{
\code{getlambda} is only applied for square matrices in which the set of \emph{n} actors is also the set of \emph{n} partners. 
}

\value{
\item{lambda}{ An \emph{n} x \emph{n} matrix of unweighted contributions to symmetry is provided. By means of this matrix dyadic and individual contributions weighted by interaction frequency will be computed.}
}

\references{

Solanas, A., Salafranca, Ll., Riba, C., Sierra, V., & Leiva, D. (2006). Quantifying social asymmetric structures. \emph{Behavior Research Methods}, \emph{38}, 390-399. 
}

\author{
David Leiva \email{dleivaur@ub.edu}.
}

\seealso{ \code{\link{getlambdaj}},\code{\link{getratiolambda}},\code{\link{reciptest1}}
} 

\examples{
##
  X=matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,ncol=3)
}

\keyword{misc}

\keyword{htest}
