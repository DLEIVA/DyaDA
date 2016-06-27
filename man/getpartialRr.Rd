\name{getpartialRr}

\alias{getpartialRr}

\title{Partial rowwise Dietz's R statistic}

\description{
Function to compute partial rowwise Dietz's R statistic.}

\usage{
   getpartialRr(X,Y,Z,names=NULL,perm=9999)
}

\arguments{
  \item{X}{Original matrix X.}
  \item{Y}{Original matrix Y.}
  \item{Z}{Original matrix Z.}
  \item{names}{Logical flag. If TRUE, a vector of individuals' names should be provided.}
  \item{perm}{Number of permutations.}
}

\details{
\code{getpartialRr} computes partial rowwise Dietz's R statistic and provides statistical significance by means of a permutation test. The function also provides several measures of the individual contributions to partial correlation. 
}

\value{

\item{rhorwxy.z}{Partial rowwise Dietz's statistic.}
\item{p_value_right}{Right-tailed p-value.}
\item{p_value_left}{Left-tailed p-value.}
\item{rhorwxy}{Correlation between X and Y partialling out effect of Z.}
\item{rhorwxz}{Correlation between X and Z partialling out effect of Y.}
\item{rhorwyz}{Correlation between Y and Z partialling out effect of X.}
\item{partial_rowwise}{Row contributions to partial rowwise Dietz's statistic}
}

\references{

de Vries, H. (1993). The rowwise correlation between two proximity matrices and the partial rowwise correlation. \emph{Psychometrika}, \emph{58}, 53-69.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.

}

\seealso{ \code{\link{getRr}}.
}

\keyword{misc}

\keyword{htest}

