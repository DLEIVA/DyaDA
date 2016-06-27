\name{getwinloss}

\alias{getwinloss}

\title{Win-loss measures at individual level}

\description{
Function to obtain win and loss measures at individual level from the observed sociomatrix.}

\usage{
   getwl(X, names=NULL, method=c("Dij","Pij"))
}

\arguments{
  \item{X}{Empirical sociomatrix containing wins-losses frequencies in dyadic encounters. The matrix must be square and numeric.}
  \item{names}{Character vector with the names of individuals. This vector is NULL by default}
  \item{method}{A character string indicating which dyadic dominance measure is to be used for the computation of David's scores. One of "Dij" or "Pij", can be abbreviated.}
}

\details{
By means of the empirical sociomatrix of wins and losses this function computes several win-loss measures at individual level. Specifically, it computes \emph{w}, \emph{weighted.w}, \emph{l} and \emph{weighted.l}. \emph{w} is the sum of individuals' dyadic dominances Dij or the sum of proportions of wins Pij by rows, depending on the specification of the \emph{method}. \emph{weighted.w} measures is the sum of individuals' Dij or Pij values weighted by the \emph{w} values of their interactants. \emph{l} is the sum of individuals' dyadic dominance indices Dij or the sum of individuals' proportions of wins Pij by columns. And finally, \emph{weighted.l} is the columns sum of individuals' Dij or Pij values weighted by the \emph{l} values of their interactants. These measures are used when computing David's scores.

}

\value{

  The result is a data frame with the following components: 
  \item{w}{Sum of dyadic dominance indices Dij or proportions of wins Pij by rows.}
  \item{weighted.w}{Weighted sum of dyadic dominance indices Dij or proportions of wins Pij.}
  \item{l}{Sum of dyadic dominance indices Dij or proportions of wins Pij by columns.}
  \item{weighted.l}{Weighted sum of dyadic dominance indices Dij or proportions of wins Pij.}
}

\references{
David, H. A. (1988). \emph{The Method of Paired Comparisons}. London: C. Griffin.

de Vries, H., Stevens, J. M. G., & Vervaecke, H. (2006). Measuring and testing the steepness of dominance hierarchies. \emph{Animal Behaviour}, \emph{71}, 585-592.
}

\author{
David Leiva \email{dleivaur@ub.edu} & Han de Vries \email{J.deVries1@uu.nl}.
}

\seealso{ \code{\link{getDS}}.
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

res <- getwinloss(X,names=individuals,method="Dij")

print(res,digits=3)

}

\keyword{misc}

