\name{getkappa}

\alias{getkappa}

\title{Dyadic reciprocity index - kappa -}

\description{
Function to obtain the dyadic reciprocity index by means of the observed sociomatrix.
}


\usage{
   getkappa(X)
}

\arguments{
  \item{X}{Original sociomatrix}
}

\details{
\code{getkappa} is only applied for square matrices in which the set of \emph{n} actors is also the set of \emph{n} partners. 
}

\value{
\item{kappa}{ The dyadic reciprocity index is computed by means of the dyadic and individual contributions to the overall skew-symmetry. It ranges from 0 to 1.}
}

\references{

Solanas, A., Salafranca, Ll., Riba, C., Sierra, V., & Leiva, D. (2006). Quantifying social asymmetric structures. \emph{Behavior Research Methods}, \emph{38}, 390-399. 
}

\author{
David Leiva \email{dleivaur@ub.edu}.
}

\seealso{ \code{\link{getnuj}}, \code{\link{getrationu}}, \code{\link{reciptest}}
} 

\examples{
##
  X=matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,ncol=3)
}

\keyword{misc}

\keyword{htest}
