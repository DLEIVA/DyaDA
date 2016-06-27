\name{getZ}

\alias{getZ}

\title{Mantel's Z statistic}

\description{
Function to compute Mantel's Z statistic.}

\usage{
   getZ(X,Y)
}

\arguments{
  \item{X}{Original matrix X.}
  \item{Y}{Original matrix Y.}
}

\details{
\code{getZ} computes Mantel's Z statistic, a measure for quantifying association between two matrices.
}

\value{

\item{Z}{Mantel's statistic.}
}

\references{

Mantel, N. (1967). The detection of disease clustering and a generalized regression approach. \emph{Cancer Research}, \emph{27}, 209-220.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.

}

\seealso{ \code{\link{getexpecZ}}, \code{\link{getSEZ}}, \code{\link{mantelZtest}}.
}

\keyword{misc}

\keyword{htest}

