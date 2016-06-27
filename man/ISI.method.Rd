\name{ISI.method}

\alias{ISI.method}

\title{I&SI method in R}

\description{
An R implementation of the I&SI algorithm for obtaining an }

\usage{
   ISI.method(X,names=NULL,tries=100)
}

\arguments{
  \item{X}{Original sociomatrix.}
  \item{names}{If TRUE, character vector with individuals' names should be given.}
  \item{tries}{Number of iterations for the I&SI algorithm.}

}

\details{
\code{ISI.method} obtains a best fitted order according to a linear or near-linear hieararchy (de Vries, 1998).
}

\references{

de Vries, H. (1998). Finding a dominance order most consistent with a linear hierarchy: a new procedure and review. \emph{Animal Behaviour}, \emph{55}, 827-843.
}

\author{
Han de Vries <J.deVries1@uu.nl>, & David Leiva <dleivaur@ub.edu>.
}

\examples{

  X <- matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,byrow=TRUE)
}

\keyword{misc}

\keyword{htest}

