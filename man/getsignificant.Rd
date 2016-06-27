\name{gettied}

\alias{gettied}

\title{Dyads with significant relationships}

\description{
Function to obtain the number of significant relationships of the observed sociomatrix.}

\usage{
   getsignificant(X,alfa=0.05,relfreq=TRUE)
}

\arguments{
  \item{X}{Original sociomatrix.}
  \item{alfa}{Significance level for the binomial test.}
  \item{relfreq}{Logical. If TRUE the function computes percentage of significant relationships.}

}

\details{
\code{gettied} computes the number of significant dyads in a binomial test at \code{alfa} level. The function carries out binomial test for every single cell of the sociomatrix using the frequency in each particular cell as the number of occurrences and the dyadic total of interactions as the total of events. 
}

\value{
\item{binom.test.matrix}{Matrix with the results of the dyadic binomial tests. Dyadic total number of interactions are shown in the upper triangular part of the matrix, whereas p-values are place in the lower triangular matrix.}
\item{tied.dyads}{Number of dyads with significant outcomes according to a binomial test.}
\item{perc.significant.dyads}{Percentage of significant dyads according to a binomial test.}
}

\references{

van Hooff, J. A. R. A. M., & Wensing, J. A. B. (1987). Dominance and its behavioral measures in a captive wolf pack. In H. W. Frank (Ed.), \emph{Man and wolf} (pp. 219-252). Dordrecht: Junk Publishers.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.
}

\seealso{ \code{\link{getonewayrel}}, \code{\link{gettied}}, \code{\link{gettwowayrel}}, \code{\link{getunknown}}.
}

\examples{

  set.seed(123)
  X <- matrix(floor(runif(20*20,min=0,max=20)),nrow=20,byrow=TRUE)
  diag(X)<-0
  result <- getsignificant(X)
  plot(result)
}

\keyword{misc}

\keyword{htest}

