\name{rowwiseRrtest}

\alias{rowwiseRrtest}

\title{Statistical significance for Rr statistic}

\description{
Function to compute statistical significance for rowwise Dietz's R statistic.}

\usage{
   rowwiseRrtest(X,Y,perm=9999)
}

\arguments{
  \item{X}{Original matrix X.}
  \item{Y}{Original matrix Y.}
  \item{perm}{Number of permutations. Between 1 and 1000000}
}

\details{
\code{rowwiseRrtest} estimates statistical significance for rowwise Dietz's R statistic by means of a permutation test. 
}

\value{

\item{Rr}{Rowwise Dietz's statistic.}
\item{right_pvalue}{Right-tailed p-value.}
\item{left_pvalue}{Left-tailed p-value}
}

\references{

de Vries, H. (1993). The rowwise correlation between two proximity matrices and the partial rowwise correlation. \emph{Psychometrika}, \emph{58}, 53-69.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.

}

\seealso{ \code{\link{getpartialRr}}, \code{\link{getRr}}.
}

\keyword{misc}

\keyword{htest}

