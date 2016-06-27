\name{block.SRM.effects}

\alias{block.SRM.effects}

\title{Estimating SRM effects in block designs}

\description{
A function for obtaining SRM effects in block designs.
}

\usage{
   block.SRM.effects(X,subgroup1,subgroup2,numb.times,numb.groups,names1=NULL,names2=NULL,times=NULL,groups=NULL)
}

\arguments{
  \item{X}{Original block sociomatrix.}
  \item{subgroup1}{Number of individuals of the first subgroup of the block data.}
  \item{subgroup2}{Number of individuals of the second subgroup of the block data.}
  \item{numb.times}{Number of occasions in which the dyads are analyzed.}
  \item{numb.groups}{Number of groups in the block data.}
  \item{names1}{ (optional) Character array with the names of individuals in the first subgroup.}
  \item{names2}{ (optional) Character array with the names of individuals in the second subgroup.}
  \item{times}{ (optional) Character array with the names of the occasions of data collection.}
  \item{groups}{ (optional) Character array with the names of groups in the block data.}
}

\details{
\code{block.SRM.effects} estimates actor,partner and relationship effects for SRM block designs.
}

\value{
\code{block.SRM.effects} returns a list containing the following components:

\item{actor.effects}{ Estimation of the actor effects for the SRM block data.}
\item{partner.effects}{ Estimation of the partner effects for the SRM block data.}
\item{relationship.effects}{ Estimation of the relationship effects for the SRM block data.}
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
