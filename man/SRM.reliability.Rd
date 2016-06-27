\name{SRM.reliability}

\alias{SRM.reliability}

\title{Estimating reliability for the SRM estimates in round robin designs}

\description{
A function to compute reliability for the estimates of the actor and partner variances in round robin designs.
}

\usage{
   SRM.reliability(X, numb.individuals, numb.times, numb.groups, names=NULL, times=NULL, groups=NULL)
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
\code{SRM.reliability} estimates reliability for the actor and partner variances in SRM round robin designs.
}

\value{
\code{SRM.reliability} returns a list containing the following components:

\item{actor.reliability}{ Estimation of reliability of the actor variance for the SRM round robin data.}
\item{partner.reliability}{ Estimation of reliability of the partner variance for the SRM round robin data.}
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
