\name{gett}

\alias{gett}

\title{Student's t statistic for two matrices}

\description{
Function to compute Student's t statistic for two matrices.}

\usage{
   gett(X,Y,method)
}

\arguments{
  \item{X}{Original matrix X.}
  \item{Y}{Original matrix Y.}
  \item{method}{String. Two values: \emph{pearson} and \emph{spearman}}
}

\details{
\code{gett} computes t statistic by means of a matrix correlation index (Pearson's or Spearman's correlation).
}

\value{

\item{t}{Student's statistic.}
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.

}

\seealso{ \code{\link{getpearson}}, \code{\link{getspearman}}.
}

\keyword{misc}

\keyword{htest}

