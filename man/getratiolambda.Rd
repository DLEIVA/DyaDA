\name{getratiolambda}

\alias{getratiolambda}

\title{Matrix of dyadic contributions to symmetry - ratiolambda -}

\description{
Function to obtain the matrix of weighted dyadic contributions to symmetry by means of the observed sociomatrix.
}


\usage{
   getratiolambda(X)
}

\arguments{
  \item{X}{Original sociomatrix}
}

\details{
\code{getratiolambda} is only applied for square matrices in which the set of \emph{n} actors is also the set of \emph{n} partners. 
}

\value{
\item{ratiolambda}{ An \emph{n} x \emph{n} matrix of dyadic weighted contributions to symmetry is provided.}
}

\references{

Solanas, A., Salafranca, Ll., Riba, C., Sierra, V., & Leiva, D. (2006). Quantifying social asymmetric structures. \emph{Behavior Research Methods}, \emph{38}, 390-399. 
}

\author{
David Leiva \email{dleivaur@ub.edu}.
}

\seealso{ \code{\link{getlambda}},\code{\link{getlambdaj}},\code{\link{reciptest1}}
} 

\examples{
##
  X=matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,ncol=3)
}

\keyword{misc}

\keyword{htest}
