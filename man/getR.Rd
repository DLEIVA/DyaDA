\name{getR}

\alias{getR}

\title{Dietz's R statistic}

\description{
Function to compute Dietz's R statistic.}

\usage{
   getR(X,Y)
}

\arguments{
  \item{X}{Original matrix X.}
  \item{Y}{Original matrix Y.}
}

\details{
\code{getR} computes Dietz's R statistic, a measure for quantifying association between two matrices.
}

\value{

\item{R}{Dietz's statistic.}
}

\references{

Dietz, E. J. (1983). Permutation tests for association between two distance matrices. \emph{Systematic Zoology}, \emph{32}, 21-26.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.

}

\seealso{ \code{\link{getexpecR}}, \code{\link{getSER}}, \code{\link{dietzRtest}}.
}

\keyword{misc}

\keyword{htest}


