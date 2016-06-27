\name{block.SRM.variances}

\alias{block.SRM.variances}

\title{Estimating SRM absolute variances in block designs}

\description{
A function for obtaining SRM absolute variances in block designs.
}

\usage{
   block.SRM.variances(X, subgroup1, subgroup2, numb.times, numb.groups, names1=NULL, names2=NULL, subgroups=NULL, times=NULL, groups=NULL)
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
\code{block.SRM.variances} estimates actor,partner and relationship variances, as well as actor-partner and relationship covariances, for SRM block designs.
}

\value{
\code{block.SRM.variances} returns a list containing the following components:

\item{actor.variance}{ Estimation of the actor variance for the SRM block data.}
\item{partner.variance}{ Estimation of the partner variance for the SRM block data.}
\item{relationship.variance}{ Estimation of the relationship variance for the SRM block data.}
\item{actorpartner.covariance}{ Estimation of the actor-partner covariance for the SRM block data.}
\item{dyadic.covariance}{ Estimation of the relationship covariance for the SRM block data.}
}

\references{

Kenny, D. A., Kashy, D. A., & Cook, W. L. (2006). \emph{Dyadic data analysis}. New York: Guilford Press.

}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, & David A. Kenny <david.kenny@uconn.edu>.
}

\seealso{ \code{\link{block.SRM.effects}},\code{\link{block.SRM.relative.variances}}
}

\keyword{htest}
