\name{block.SRM.dyadic.reciprocity}

\alias{block.SRM.dyadic.reciprocity}

\title{Computing dyadic reciprocity in SRM block designs}

\description{
A function for computing dyadic reciprocity in SRM block designs.
}


\usage{
   block.SRM.dyadic.reciprocity(X,subgroup1,subgroup2,numb.times,numb.groups,symmetric=FALSE,names1=NULL,names2=NULL,subgroups=NULL,times=NULL,groups=NULL)
}

\arguments{
  \item{X}{Original block sociomatrix.}
  \item{subgroup1}{Number of individuals of the first subgroup of the block data.}
  \item{subgroup2}{Number of individuals of the second subgroup of the block data.}
  \item{numb.times}{Number of occasions in which the dyads are analyzed.}
  \item{numb.groups}{Number of groups in the block data.}
  \item{symmetric}{If 0 asymmetrical block designs are assumed. Otherwise, formulas for symmetrical block designs are used.}
  \item{names1}{ (optional) Character array with the names of individuals in the first subgroup.}
  \item{names2}{ (optional) Character array with the names of individuals in the second subgroup.}
  \item{subgroups}{ (optional) Character array with the names of the subgroups.}
  \item{times}{ (optional) Character array with the names of the occasions of data collection.}
  \item{groups}{ (optional) Character array with the names of groups in the block data.}
}

\details{
\code{block.SRM.dyadic.reciprocity} estimates dyadic reciprocity for symmetric and non-symmetric SRM block designs.
}

\value{
\item{dyadic.reciprocity}{Dyadic reciprocity measure for SRM block designs.}
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
