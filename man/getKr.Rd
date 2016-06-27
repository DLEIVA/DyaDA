\name{getKr}

\alias{getKr}

\title{Rowwise Kendall's K statistic -Kr-}

\description{
Function to compute rowwise Kendall's K statistic.}

\usage{
   getKr(X,Y,names=NULL)
}

\arguments{
  \item{X}{Original matrix X.}
  \item{Y}{Original matrix Y.}
  \item{names}{Logical flag. If TRUE, a vector of individuals' names should be provided.}
}

\details{
\code{getKr} computes rowwise Kendall's K statistic. The function also provides several measures of the individual contributions to Kendall's tau and Kr statistics. 
}

\value{

\item{Kr}{Rowwise Kendall's statistic.}
\item{rtauwav}{Average rowwise Kendall's tau statistic.}
\item{rtauw}{Rowwise Kendall's tau statistic.}
\item{rowtau}{Row contributions to rowwise Kendall's tau statistic.}
\item{weighted_contributions}{Weighted contributions to rowwise Kendall's statistic}
\item{contributions}{Row contributions to rowwise Kendall's statistic}
}

\references{

de Vries, H. (1993). The rowwise correlation between two proximity matrices and the partial rowwise correlation. \emph{Psychometrika}, \emph{58}, 53-69.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.

}

\seealso{ \code{\link{getpartialKr}}, \code{\link{rowwiseKrtest}}.
}

\keyword{misc}

\keyword{htest}

