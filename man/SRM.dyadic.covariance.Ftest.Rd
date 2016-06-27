\name{SRM.dyadic.covariance.Ftest}

\alias{SRM.dyadic.covariance.Ftest}

\title{Testing dyadic covariance in SRM round robin designs by means of a F test}

\description{
A function for testing dyadic reciprocity in SRM round robin designs.
}

\usage{
   SRM.dyadic.covariance.Ftest(X, numb.individuals, numb.times, numb.groups, names=NULL, times=NULL, groups=NULL)
}

\arguments{
  \item{X}{Original round robin sociomatrix.}
  \item{numb.individuals}{Number of individuals in the round robin data.}
  \item{numb.times}{Number of occasions in which the dyads are analyzed.}
  \item{numb.groups}{Number of groups in the round robin data.}
  \item{names}{ (optional) Character array with the names of individuals in the group.}
  \item{times}{ (optional) Character array with the names of the occasions of data collection.}
  \item{groups}{ (optional) Character array with the names of groups in the round robin data.}
}

\details{
\code{SRM.dyadic.covariance.Ftest} estimates dyadic reciprocity for SRM round robin designs by means of a ratio of variances F test.
}

\value{
\code{SRM.dyadic.covariance.Ftest} returns a list containing the following components:

\item{F.value}{ Estimation of the ratio of variances for testing dyadic reciprocity for round robin designs.}
\item{df1}{ Degrees of freedom of the numerator.}
\item{df2}{ Degrees of freedom of the numerator.}
\item{p.value}{ Statistical significance of the empirical F statistic for the SRM round robin data.}
}

\references{

Kenny, D. A., Kashy, D. A., & Cook, W. L. (2006). \emph{Dyadic data analysis}. New York: Guilford Press.

}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, & David A. Kenny <david.kenny@uconn.edu>.
}

\seealso{ \code{\link{SRM.variances}}
}

\keyword{htest}
