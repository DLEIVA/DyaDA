\name{indistinguishable.dyad}

\alias{indistinguishable.dyad}

\title{Estimating nonindependence in indistinguishable standard dyadic designs with interval outcomes}

\description{
A function for estimating dyadic interdependence in standard dyadic designs with indistinguishable dyad members and interval responses.
}

\usage{
   indistinguishable.dyad(dataset, conf.int=0.95)
}

\arguments{
  \item{dataset}{ Original data with dyadic responses.}
  \item{conf.int}{ Confidence level in order to estimate the confidence interval for the rho parameter. It is equal to 0.95 by default}
}

\details{
\code{indistinguishable.dyad} estimates dyadic interdependence for any standard dyadic design in which indistinguishable members and interval data are taken into account.
}

\value{
\code{indistinguishable.dyad} returns a list containing the following components:

\item{descriptive.statistics}{ Summary of some descriptive statistics for the dyadic variables.}
\item{intraclass.correlation}{ ICC coefficient for the original data.}
\item{between.dyads.MS}{ Between dyads mean of sum of squares.}
\item{within.dyads.MS}{ Within dyads mean of sum of squares.}
\item{F.statistic}{ Statistic value for carrying out the statistical test regarding dyadic nonindependence.}
\item{df1}{ Degrees of freedom of the numerator.}
\item{df2}{ Degrees of freedom of the denominator.}
\item{two.tailed.p.value}{ Statistical significance of the t statistic value under the null hypothesis of no dyadic correlation}
\item{alpha}{ Probability of the parameter not included in the interval estimation.}
\item{ICC.confidence.interval}{ Lower and upper bounds of the ICC parameter estimated at 1-alpha confidence level.}
}

\references{

Kenny, D. A., Kashy, D. A., & Cook, W. L. (2006). \emph{Dyadic data analysis}. New York: Guilford Press.

}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, & David A. Kenny <david.kenny@uconn.edu>.
}

\seealso{ \code{\link{distinguishable.dyad}}
}

\keyword{htest}

