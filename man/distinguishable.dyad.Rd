\name{distinguishable.dyad}

\alias{distinguishable.dyad}

\title{Estimating nonindependence in distinguishable standard dyadic designs with interval outcomes}

\description{
A function for estimating dyadic interdependence in standard dyadic designs with distinguishable dyad members and interval responses.
}

\usage{
   distinguishable.dyad(dataset, conf.int=0.95)
}

\arguments{
  \item{dataset}{ Original data with dyadic responses.}
  \item{conf.int}{ Confidence level in order to estimate the confidence interval for the rho parameter. It is equal to 0.95 by default}
}

\details{
\code{distinguishable.dyad} estimates dyadic interdependence for any standard dyadic design in which distinguishable members and interval data are taken into account.
}

\value{
\code{distinguishable.dyad} returns a list containing the following components:

\item{descriptive.statistics}{ Summary of some descriptive statistics for the dyadic variables.}
\item{pearson.correlation}{ Pearson's correlation coefficient for the original data.}
\item{alpha}{ Probability of the parameter not included in the interval estimation.}
\item{r.pearson.confidence.interval}{ Lower and upper bounds of the rho parameter estimated at 1-alpha confidence level.}
\item{t.statistic}{ Statistic value for carrying out the statistical test regarding dyadic nonindependence.}
\item{df}{ Degrees of freedom.}
\item{two.tailed.p.value}{ Statistical significance of the t statistic value under the null hypothesis of no dyadic correlation}

}

\references{

Kenny, D. A., Kashy, D. A., & Cook, W. L. (2006). \emph{Dyadic data analysis}. New York: Guilford Press.

}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, & David A. Kenny <david.kenny@uconn.edu>.
}

\seealso{ \code{\link{indistinguishable.dyad}}
}

\keyword{htest}


