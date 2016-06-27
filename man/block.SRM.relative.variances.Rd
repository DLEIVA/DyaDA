\name{block.SRM.relative.variances}

\alias{block.SRM.relative.variances}

\title{Estimating SRM relative variances in block designs}

\description{
A function for obtaining SRM relative variances in block designs.
}

\usage{
   block.SRM.relative.variances(X, subgroup1, subgroup2, numb.times, numb.groups, names1=NULL, names2=NULL, subgroups=NULL, times=NULL, groups=NULL)
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
\code{block.SRM.relative.variances} estimates relative actor,partner and relationship variances for SRM block designs.
}

\value{
\code{block.SRM.relative.variances} returns a list containing the following components:

\item{relative.actor.variance}{ Estimation of the relative actor variance for the SRM block data.}
\item{relative.partner.variance}{ Estimation of the relative partner variance for the SRM block data.}
\item{relative.relationship.variance}{ Estimation of the relative relationship variance for the SRM block data.}
}

\references{

Kenny, D. A., Kashy, D. A., & Cook, W. L. (2006). \emph{Dyadic data analysis}. New York: Guilford Press.

}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, & David A. Kenny <david.kenny@uconn.edu>.
}

\seealso{ \code{\link{block.SRM.effects}},\code{\link{block.SRM.variances}}
}

\keyword{htest}

