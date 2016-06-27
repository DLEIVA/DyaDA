\name{SRM.generalized.reciprocity}

\alias{SRM.generalized.reciprocity}

\title{Computing generalized reciprocity in SRM round robin designs}

\description{
A function for computing generalized reciprocity in SRM round robin designs.
}

\usage{
   SRM.generalized.reciprocity(X, numb.individuals, numb.times, numb.groups, names=NULL, times=NULL, groups=NULL)
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
\code{SRM.generalized.reciprocity} estimates generalized reciprocity for SRM round robin designs.
}

\value{
\item{generalized.reciprocity}{Generalized reciprocity measure for SRM round robin designs.}
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
