\name{dietzRtest}

\alias{dietzRtest}

\title{Statistical test for Dietz's R statistic}

\description{
Function to estimate statistical significance for Dietz's R statistic.}

\usage{
   dietzRtest(X,Y,perm=9999)
}

\arguments{
  \item{X}{Original matrix X.}
  \item{Y}{Original matrix Y.}
  \item{perm}{Number of permutations. Between 1 and 1000000}
}

\details{
\code{dietzRtest} estimates statistical significance of Dietz's R statistic by means of a permutation test.
}

\value{

\item{Dietz_R}{Dietz's statistic.}
\item{right_pvalue}{Right-tailed p-value.}
\item{left_pvalue}{Left-tailed p-value}
}

\references{

Dietz, E. J. (1983). Permutation tests for association between two distance matrices. \emph{Systematic Zoology}, \emph{32}, 21-26.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.

}

\seealso{ \code{\link{getexpecR}}, \code{\link{getR}}, \code{\link{getSER}}.
}

\keyword{misc}

\keyword{htest}

