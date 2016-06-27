\name{getDS}

\alias{getDS}

\title{David's scores -DS-}

\description{
Function to obtain David's scores from the observed sociomatrix.
}

\usage{
   getDS(X, names=NULL, method=c("Dij","Pij"))
}

\arguments{
  \item{X}{Empirical sociomatrix containing wins-losses frequencies in dyadic encounters. The matrix must be square and numeric.}
  \item{names}{Character vector with the names of individuals. This vector is NULL by default}
  \item{method}{A character string indicating which dyadic dominance measure is to be used for the computation of David's scores. One of "Dij" or "Pij", can be abbreviated.}
}

\details{
\code{getDS} is obtained by means of the following expression: \eqn{DS = w1 + w2 - l1 - l2} where \emph{w1} is the sum of i's Dij or Pij values (depending on the \code{method} specification); \emph{w2} is the weighted sum of i`s dyadic dominance indices corrected for chance or the weighted sum of i's win proportions; \emph{l1} is the sum of i's Dji or Pji values and \emph{l2} is the sum of i's dyadic lose indices corrected for chance or the weighted sum of i's lose proportions.
}

\value{

\item{DS}{David's scores based on dyadic dominance indices corrected for chance or on win proportions.}
}

\references{
David, H. A. (1988). \emph{The Method of Paired Comparisons}. London: C. Griffin.

de Vries, H., Stevens, J. M. G., & Vervaecke, H. (2006). Measuring and testing the steepness of dominance hierarchies. \emph{Animal Behaviour}, \emph{71}, 585-592.
}

\author{
David Leiva \email{dleivaur@ub.edu} & Han de Vries \email{J.deVries1@uu.nl}.
}

\seealso{\code{\link{getwinloss}}.
}

\examples{

##############################################################################
###               Example taken from Vervaecke et al. (2007):              ###
##############################################################################

X <- matrix(c(0,58,50,61,32,37,29,39,25,8,0,22,22,9,27,20,10,48,
              3,3,0,19,29,12,13,19,8,5,8,9,0,33,38,35,32,57,
              4,7,9,1,0,28,26,16,23,4,3,0,0,6,0,7,6,12,
              2,0,4,1,4,4,0,5,3,0,2,1,1,5,8,3,0,10,3,1,3,0,0,4,1,2,0),
              nrow=9,byrow=TRUE)

individuals <- c("V","VS","B","FJ","PR","VB","TOR","MU","ZV")

res <- getDS(X,names=individuals,method="Dij")

print(res,digits=3)
}

\keyword{misc}

\keyword{htest}

