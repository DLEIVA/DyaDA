\name{getRr}

\alias{getRr}

\title{Rowwise Dietz's R statistic -Rr-}

\description{
Function to compute rowwise Dietz's R statistic.}

\usage{
   getRr(X,Y,names=NULL)
}

\arguments{
  \item{X}{Original matrix X.}
  \item{Y}{Original matrix Y.}
  \item{names}{Logical flag. If TRUE, a vector of individuals' names should be provided.}
}

\details{
\code{getRr} computes rowwise Dietz's R statistic. The function also provides several measures of the individual contributions to Spearman's and Dietz's statistics. 
}

\value{

\item{Rr}{Rowwise Dietz's statistic.}
\item{rrhowav}{Average rowwise Spearman's correlation.}
\item{rrhow}{Rowwise Spearman's correlation.}
\item{rowrho}{Row contributions to rowwise Spearman's statistic.}
\item{weighted_contributions}{Weighted contributions to rowwise Dietz's statistic}
\item{contributions}{Row contributions to rowwise Dietz's statistic}
}

\references{

de Vries, H. (1993). The rowwise correlation between two proximity matrices and the partial rowwise correlation. \emph{Psychometrika}, \emph{58}, 53-69.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.

}

\seealso{ \code{\link{getpartialRr}}, \code{\link{rowwiseRrtest}}.
}

\keyword{misc}

\keyword{htest}

