\name{SRM.between.groups.tTest}

\alias{SRM.between.groups.tTest}

\title{Testing SRM effects in round robin designs by means of a between-group t Test}

\description{
A function for testing SRM effects in round robin designs by means of a between-group t Test.
}


\usage{
   SRM.between.groups.tTest(X, numb.individuals, numb.times, numb.groups, names=NULL, times=NULL, groups=NULL)
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
\code{SRM.between.groups.tTest} estimates statistical significance for the variance and covariance components obtained for a SRM round robin design.
}

\value{
\code{SRM.between.groups.tTest} returns a list containing the following components:

\item{t.statistic}{ t statistic for each variance component of the SRM round robin data.}
\item{df}{Degrees of freedom for the statistical test.}
\item{two.tailed.p.value}{Statistical significance for the variance components of the SRM round robin data.}
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
