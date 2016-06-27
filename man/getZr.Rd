\name{getZr}

\alias{getZr}

\title{Rowwise Mantel's Z statistic -Zr-}

\description{
Function to compute rowwise Mantel's Z statistic.}

\usage{
   getZr(X,Y,names=NULL)
}

\arguments{
  \item{X}{Original matrix X.}
  \item{Y}{Original matrix Y.}
  \item{names}{Logical flag. If TRUE, a vector of individuals' names should be provided.}
}

\details{
\code{getZr} computes rowwise Mantel's Z statistic. The function also provides several measures of the individual contributions to Pearson's and Mantel's statistics. 
}

\value{

\item{Zr}{Rowwise Mantel's statistic.}
\item{rrwav}{Average rowwise Pearson's correlation.}
\item{rrw}{Rowwise Pearson's correlation.}
\item{rowpearson}{Row contributions to rowwise Pearson's statistic.}
\item{weighted_contributions}{Weighted contributions to rowwise Mantel's statistic}
\item{contributions}{Row contributions to rowwise Mantel's statistic}
}

\references{

de Vries, H. (1993). The rowwise correlation between two proximity matrices and the partial rowwise correlation. \emph{Psychometrika}, \emph{58}, 53-69.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.

}

\seealso{ \code{\link{getpartialZr}}, \code{\link{rowwiseZrtest}}.
}

\keyword{misc}

\keyword{htest}

