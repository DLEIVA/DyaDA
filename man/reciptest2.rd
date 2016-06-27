\name{reciptest2}

\alias{reciptest2}

\title{Statistical tests for social reciprocity}

\description{
Estimating statistical significance for social reciprocity statistics at different levels of analysis.}

\usage{
   reciptest(X, pi ,rep=9999, names=NULL, label=FALSE)
}

\arguments{
  \item{X}{Original sociomatrix. Matrix X must be square.}
  \item{pi}{Matrix of probabilities \emph{Pij} that specifies the null hypothesis regarding social reciprocity. Pi must be a square matrix}
  \item{rep}{Number of simulations for carrying out the randomization test. Number of simulations must be between 1 and 1000000.}
  \item{names}{Character vector with the names of individuals. IF \emph{label} TRUE. then this vector should be specified.}
  \item{label}{Logical flag. TRUE. if \emph{names} of individuals will be provided. Otherwise FALSE.}
}

\details{
\code{reciptest} estimates statistical significance for several measures of social reciprocity at different levels of analysis: overall measures of social reciprocity as \emph{PHIr}, see \code{\link{getPHIr}}; dyadic measures as \emph{PHIij}, see \code{\link{getPHIij}}; and individual measures as \emph{phii}, see \code{\link{getphii}}. This procedure simulates a number of sociomatrices under specified null hypothesis regarding social reciprocity (e.g. complete reciprocation) by means of callings to a C routine included in the package, then computes the reciprocity measures at different levels. After \emph{rep} simulations the sampling distributions for the statistics are estimated. Then statistical significance is computed as follows:
\eqn{p=NS+1/NOS+1}
where \emph{NS} is the number of times that simulated values is as great as or greater than the empirical values and \emph{NOS} represents the number of simulated values when \emph{p-right} value is computed. When estimating \emph{p-left} value \emph{NS} is the number of times that simulated values is as great as or lower than the empirical values.
}

\value{

\item{PHIr}{Overall measure of asymmetry in social interactions.}
\item{PHIrexpec}{Mathematical expectancy of the PHIr statistic.}
\item{PHIrSE}{Standard error of the PHIr statistic.}
\item{phii}{Individuals' contributions to asymmetry as actors.}
\item{phij}{Individuals' contributions to asymmetry as partners.}
\item{phirmat}{Dyadic directional contributions to asymmetry.}
\item{PHIij}{Dyadic contributions to asymmetry.}
}

\references{

Solanas, A., Leiva, D., Sierra, V., & Salafranca, Ll. (2009). Measuring and making decisions for social reciprocity. \emph{Behavior Research Methods}, \emph{41}, 742-754.
}

\author{
David Leiva \email{dleivaur@ub.edu}, Antonio Solanas \email{antonio.solanas@ub.edu}.
}

\seealso{ \code{\link{getPHIr}}, \code{\link{getphii}}, \code{\link{getphij}}, \code{\link{getphirmat}}, \code{\link{getPHIij}}
}

\examples{
  X=matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,ncol=3)
  pi=matrix(c(0, 0.5, 0.5, 0.5, 0, 0.5, 0.5, 0.5, 0),nrow=3,ncol=3)
  rep=99999
  names=c("Ind.1","Ind.2","Ind.3")  
}

\keyword{misc}

\keyword{htest}
