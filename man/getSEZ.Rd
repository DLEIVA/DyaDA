\name{getSEZ}

\alias{getSEZ}

\title{Standard Error for Mantel's Z statistic}

\description{
Function to compute standard error for Mantel's Z statistic.}

\usage{
   getSEZ(X,Y)
}

\arguments{
  \item{X}{Original matrix X.}
  \item{Y}{Original matrix Y.}
}

\details{
\code{getSEZ} computes standard error for Mantel's Z statistic.
}

\value{

\item{SEZ}{Standard error for Mantel's Z statistic.}
}

\references{

Mantel, N. (1967). The detection of disease clustering and a generalized regression approach. \emph{Cancer Research}, \emph{27}, 209-220.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.

}

\seealso{ \code{\link{getexpecZ}}, \code{\link{getZ}}, \code{\link{mantelZtest}}.
}

\keyword{misc}

\keyword{htest}

