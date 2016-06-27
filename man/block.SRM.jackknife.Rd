\name{block.SRM.jackknife}

\alias{block.SRM.jackknife}

\title{Testing SRM effects in block designs by means of Jackknife tests}

\description{
A function for testing SRM effects in block designs by means of a Jackknife test.
}


\usage{
   block.SRM.jackknife(X, subgroup1, subgroup2, numb.times, numb.groups, names1=NULL, names2=NULL, subgroups=NULL, times=NULL, groups=NULL)
}

\arguments{
  \item{X}{Original block sociomatrix.}
  \item{subgroup1}{Number of individuals of the first subgroup of the block data.}
  \item{subgroup2}{Number of individuals of the second subgroup of the block data.}
  \item{numb.times}{Number of occasions in which the dyads are analyzed.}
  \item{numb.groups}{Number of groups in the block data.}
  \item{names1}{ (optional) Character array with the names of individuals in the first subgroup.}
  \item{names2}{ (optional) Character array with the names of individuals in the second subgroup.}
  \item{subgroups}{ (optional) Character array with the names of the subgroups.}
  \item{times}{ (optional) Character array with the names of the occasions of data collection.}
  \item{groups}{ (optional) Character array with the names of groups in the block data.}
}

\details{
\code{block.SRM.jackknife} estimates statistical significance for the variance and covariance components obtained for a SRM block design.
}

\value{
\code{block.SRM.jackknife} returns a list containing the following components:

\item{t.statistic}{ t statistic for each variance component of the SRM block data estimated by means of a Jackknife procedure.}
\item{df}{Degrees of freedom for the statistical test.}
\item{two.tailed.p.value}{Statistical significance for the variance components of the SRM block data.}
}

\references{

Kenny, D. A., Kashy, D. A., & Cook, W. L. (2006). \emph{Dyadic data analysis}. New York: Guilford Press.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, & David A. Kenny <david.kenny@uconn.edu>.
}

\seealso{ \code{\link{block.SRM.variances}}
}

\keyword{htest}


