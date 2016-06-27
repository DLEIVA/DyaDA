\name{gettaukendall}

\alias{gettaukendall}

\title{Kendall's statistics proposed by Dietz}

\description{
Function to compute Kendall's statistics proposed by Dietz.}

\usage{
   gettaukendall(X,Y)
}

\arguments{
  \item{X}{Original matrix X.}
  \item{Y}{Original matrix Y.}
}

\details{
\code{gettaukendall} computes several decompositions of Kendall's tau statistic proposed by Dietz (1983).
}

\value{
\item{K}{Kendall's tau statistic.}
\item{Kc}{Kendall's tau statistic.}
\item{Kr}{Kendall's tau statistic.}
\item{Ku}{Kendall's tau statistic.}
}

\references{

Dietz, E. J. (1983). Permutation tests for association between two distance matrices. \emph{Systematic Zoology}, \emph{32}, 21-26.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.

}

\seealso{ \code{\link{getKr}}.
}

\keyword{misc}

\keyword{htest}


