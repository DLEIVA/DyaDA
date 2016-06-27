\name{reciptest1}
\alias{reciptest1}
\title{Social reciprocity measures and tests at different levels of analysis}
\description{
Statistical significance for social reciprocity at different levels of analysis.}
\usage{
   reciptest(X, pi, rep=9999, names=NULL, label=FALSE)
}
\arguments{
  \item{X}{Original sociomatrix. Matrix X must be square.}
  \item{pi}{Matrix of probabilities \emph{Pij} that specifies the null hypothesis regarding social reciprocity. Pi must be a square matrix}
  \item{rep}{Number of simulations for carrying out the randomization test. Number of simulations must be between 1 and 1000000.}
  \item{names}{ (optional) Character array with the names of individuals. IF \code{label} TRUE. then this array should be specified.}
  \item{label}{  (optional) Logical flag. TRUE. if \code{names} of individuals will be provided.}
}
\details{
\code{reciptest} estimates statistical significance for several measures of social reciprocity at different levels of analysis: overall measures of social reciprocity as \emph{phi}, see \code{\link{getphi}}; dyadic measures as \emph{ratiolambda}, see \code{\link{getratiolambda}}; and individual measures as \emph{nuj}, see \code{\link{getrationu}}. This procedure simulates a number of sociomatrices under specified null hypothesis regarding social reciprocity (e.g. complete reciprocation) by means of callings to a C routine included in the package, then computes the reciprocity measures at different levels. After \emph{rep} simulations the sampling distributions for the statistics are estimated. Then statistical significance is computed as follows:
\eqn{p=NS+1/NOS+1}
where \emph{NS} is the number of times that simulated values is as great as or greater than the empirical values and \emph{NOS} represents the number of simulated values when \emph{p-right} value is computed. When estimating \emph{p-left} value \emph{NS} is the number of times that simulated values is as great as or lower than the empirical values.
}
\value{
\item{phi}{Overall measure of skew-symmetry in social interactions.}
\item{psi}{Overall measure of symmetry in social interactions.}
\item{dc}{Directional consistency in social interactions.}
\item{delta}{Ratio skew-symmetry/symmetry.}
\item{epsilon}{Generalized reciprocity index in social interactions.}
\item{kappa}{Dyadic reciprocity index in social interactions.}
\item{nuj}{Individuals' contributions to skew-symmetry.}
\item{lambdaj}{Individuals' contributions to symmetry.}
\item{rationu}{Dyadic contributions to skew-symmetry.}
\item{ratiolambda}{Dyadic contributions to symmetry.}
\item{omega}{Dyadic balanced reciprocity matrix.}
}
\references{
Leiva, D., Solanas, A., & Salafranca, Ll. (2008). Testing reciprocity in social interactions: A comparison between the directional consistency and skew-symmetry statistics. \emph{Behavior Research Methods}, \emph{40}, 626-634. 

Solanas, A., Salafranca, Ll., Riba, C., Sierra, V., & Leiva, D. (2006). Quantifying social asymmetric structures. \emph{Behavior Research Methods}, \emph{38}, 390-399. 

van Hooff, J. A. R. A. M., & Wensing, J. A. B. (1987). Dominance and its behavioral measures in a captive wolf pack. In H. W. Frank (Ed.), \emph{Man and wolf} (pp. 219-252). Dordrecht: Junk Publishers.
}
\author{
David Leiva \email{dleivaur@ub.edu}.
}
\seealso{ \code{\link{getphi}}, \code{\link{getdc}}, \code{\link{getnuj}}, \code{\link{getrationu}}, \code{\link{getomega}}
} 
\examples{
  X=matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,ncol=3)
  pi=matrix(c(0, 0.5, 0.5, 0.5, 0, 0.5, 0.5, 0.5, 0),nrow=3,ncol=3)
  rep=99999
  names=c("Ind1","Ind2","Ind3")  
}
\keyword{misc}
\keyword{htest}
