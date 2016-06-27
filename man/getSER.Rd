\name{getSER}

\alias{getSER}

\title{Standard Error for Dietz's R statistic}

\description{
Function to compute standard error for Dietz's R statistic.}

\usage{
   getSER(X,Y)
}

\arguments{
  \item{X}{Original matrix X.}
  \item{Y}{Original matrix Y.}
}

\details{
\code{getSER} computes standard error for Dietz's R statistic.
}

\value{

\item{SER}{Standard error for Dietz's Z statistic.}
}

\references{

Dietz, E. J. (1983). Permutation tests for association between two distance matrices. \emph{Systematic Zoology}, \emph{32}, 21-26.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.

}

\seealso{ \code{\link{getexpecR}}, \code{\link{getR}}, \code{\link{dietzRtest}}.
}

\keyword{misc}

\keyword{htest}

