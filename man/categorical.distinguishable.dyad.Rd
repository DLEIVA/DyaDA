\name{categorical.distinguishable.dyad}

\alias{categorical.distinguishable.dyad}

\title{Estimating nonindependence in distinguishable standard dyadic designs with categorical outcomes}

\description{
A function for estimating dyadic interdependence in standard dyadic designs with distinguishable dyad members and categorical responses.
}

\usage{
   categorical.distinguishable.dyad(dataset, conf.int=0.95)
}

\arguments{
  \item{dataset}{ Original data with dyadic responses.}
  \item{conf.int}{ Confidence level in order to estimate the confidence interval for the kappa parameter. It is equal to 0.95 by default.}
}

\details{
\code{categorical.distinguishable.dyad} estimates dyadic interdependence for any standard dyadic design in which distinguishable members and categorical data are taken into account.
}

\value{
\code{categorical.distinguishable.dyad} returns a list containing the following components:

\item{contingency.table}{ Crosstable for the categorical dyadic data.}
\item{kappa}{ Cohen's kappa statistic for the original data.}
\item{alpha}{ Probability of the parameter not included in the interval estimation.}
\item{confidence.interval}{ Lower and upper bounds of the kappa parameter estimated at 1-alpha confidence level.}
\item{z.value}{ Statistic value for carrying out the statistical test regarding dyadic nonindependence.}
\item{standard.error}{ Standard error for kappa parameter under the null hypothesis of no concordance between dyad members.}
\item{two.tailed.p.value}{ Statistical significance of the z value under the null hypothesis of no concordance between dyad members}

}

\references{

Kenny, D. A., Kashy, D. A., & Cook, W. L. (2006). \emph{Dyadic data analysis}. New York: Guilford Press.

}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, & David A. Kenny <david.kenny@uconn.edu>.
}

\seealso{ \code{\link{categorical.indistinguishable.dyad}}
}

\keyword{htest}

