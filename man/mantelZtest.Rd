\name{mantelZtest}

\alias{mantelZtest}

\title{Statistical test for Mantel's Z statistic}

\description{
Function to estimate statistical significance for Mantel's Z statistic.}

\usage{
   mantelZtest(X,Y,perm=9999)
}

\arguments{
  \item{X}{Original matrix X.}
  \item{Y}{Original matrix Y.}
  \item{perm}{Number of permutations. Between 1 and 1000000}
}

\details{
\code{mantelZtest} estimates statistical significance of Mantel's Z statistic by means of a permutation test.
}

\value{

\item{mantel_Z}{Mantel's statistic.}
\item{right_pvalue}{Right-tailed p-value.}
\item{left_pvalue}{Left-tailed p-value}
}

\references{

Mantel, N. (1967). The detection of disease clustering and a generalized regression approach. \emph{Cancer Research}, \emph{27}, 209-220.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.

}

\seealso{ \code{\link{getexpecZ}}, \code{\link{getSEZ}}, \code{\link{getZ}}.
}

\keyword{misc}

\keyword{htest}

