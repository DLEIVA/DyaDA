\name{getpvkendall}

\alias{getpvkendall}

\title{Statistical Significance for Kendall's Linearity Index}

\description{
Function to estimate statistical significance for Kendall's linearity index.}

\usage{
   getpvkendall(X)
}

\arguments{
  \item{X}{Original sociomatrix.}

}

\details{
\code{getpvkendall} estimates statistical significance for Kendall's linearity index by means of a Chi-squared test. It is appropriated only for groups of size bigger than or equal to 10 individuals.
}

\value{

\item{Chisq}{Chi-squared statistical value.}
\item{df}{Degrees of freedom.}
\item{pvalue}{Statistical significance for the statistic.}
}

\references{

Kendall, M. G., & Babington Smith, B. (1940). On the method of paired comparisons. \emph{Biometrika}, \emph{31}, 324-345.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.
}

\seealso{ \code{\link{getcirc}}, \code{\link{getexpeccirc}}, \code{\link{getkendall}}, \code{\link{getmaxcirc}}.
}

\examples{

  set.seed(123)
  X <- matrix(floor(runif(20*20,min=0,max=20)),nrow=20,byrow=TRUE)
  diag(X)<-0
  getpvkendall(X)
}

\keyword{misc}

\keyword{htest}

