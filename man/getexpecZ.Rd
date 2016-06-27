\name{getexpecZ}

\alias{getexpecZ}

\title{Mathematical Expectancy for Mantel's Z statistic}

\description{
Function to compute mathematical expectancy for Mantel's Z statistic.}

\usage{
   getexpecZ(X,Y)
}

\arguments{
  \item{X}{Original matrix X.}
  \item{Y}{Original matrix Y.}
}

\details{
\code{getexpecZ} computes mathematical expectancy for Mantel's Z statistic.
}

\value{

\item{expecZ}{Mathematical Expectancy for Mantel's Z statistic.}
}

\references{

Mantel, N. (1967). The detection of disease clustering and a generalized regression approach. \emph{Cancer Research}, \emph{27}, 209-220.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.

}

\seealso{ \code{\link{getSEZ}}, \code{\link{getZ}}, \code{\link{mantelZtest}}.
}

\keyword{misc}

\keyword{htest}



