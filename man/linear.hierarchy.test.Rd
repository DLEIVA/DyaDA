\name{linear.hierarchy.test}

\alias{linear.hierarchy.test}

\title{Testing Linearity in dominanace hierarchies}

\description{
Function to estimate statistical significance for Landau's h index.}

\usage{
   linear.hierarchy.test(X,rep=9999)
}

\arguments{
  \item{X}{Original sociomatrix.}
  \item{rep}{Number of randomizations. A number between 1 and 1000000 should be specified.}

}

\details{
\code{linear.hierarchy.test} estimates statistical significance for Landau's h index by means of a two-step randomization procedure proposed by de Vries (1995). If number of individuals is greater or equal to 10 it also provides with statistical significance for Kendall's linearity index by means of a Chi-squared test.
}

\value{
\item{nind}{Number of individuals.}
\item{dyads}{Number of dyads.}
\item{total}{Total number of interactions.}
\item{desctable}{Description of the dyadic relationships in the input sociomatrix.}
\item{randomizations}{Specified number of randomizations.}
\item{landauh}{Landau's h index.}
\item{improved.landauh}{Improved Landau's index h'.}
\item{expected.h}{Mathematical expectancy of Landau's h under the assumption of no linear hierarchy.}
\item{variance.h}{Variance of Landau's h under the assumption of no linear hierarchy.}
\item{right.pvalue}{Right-tailed p-value.}
\item{left.pvalue}{Left-tailed p-value.}
\item{circ}{Number of circular triads.}
\item{expcirc}{Expected number of circular triads.}
\item{maxcirc}{Maximum number of circular triads.}
\item{rho}{Kendall's rho index.}
\item{rho.test}{Statistical test for Kendall's rho.}
}

\references{

de Vries, H. (1995). An improved test of linearity in dominance hierarchies containing unknown or tied relationships. \emph{Animal Behaviour}, \emph{50}, 1375-1389.

Kendall, M. G., & Babington Smith, B. (1940). On the method of paired comparisons. \emph{Biometrika}, \emph{31}, 324-345.

}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.
}

\seealso{ \code{\link{getexpech}}, \code{\link{getimplandau}}, \code{\link{getlandau}}, \code{\link{getvarh}}, \code{\link{getcirc}}, \code{\link{getexpeccirc}}, \code{\link{getkendall}}, \code{\link{getmaxcirc}}.
}

\examples{

  set.seed(123)
  X <- matrix(floor(runif(20*20,min=0,max=20)),nrow=20,byrow=TRUE)
  diag(X)<-0
  linear.hierarchy.test(X)
}

\keyword{misc}

\keyword{htest}

