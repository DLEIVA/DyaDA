\name{freq2abil}

\alias{freq2abil}

\title{Transformation of a matrix of social interactions}

\description{
Function to transform a sociomatrix into a matrix of dominance according.}

\usage{
   freq2abil(X,method=c("landau","matman","schmid"))
}

\arguments{
  \item{X}{Original sociomatrix.}
  \item{method}{String vector. 3 possible values: \code{"landau"}, \code{"matman"}, and \code{"schmid"}.}
}

\details{
\code{freq2abil} transforms the matrix of agonistic interactions into dyadic dominance scores according to 3 possibilities:

  - \code{"landau"}: A matrix of -1s and 1s is obtained for submissive and dominant relationships, respectively. Null and tied dyads produce 0s in the matrix of dominance scores.
 
  - \code{"matman"}: A matrix of 1s and 0s is obtained for submissive and dominant relationships, respectively. Null and tied dyads produce 0.5s in the matrix of dominance scores.

  - \code{"schmid"}: A matrix of -1s and 1s is obtained for submissive and dominant relationships, respectively. Null and tied dyads respectively produce 0s and 0.5s in the matrix of dominance scores.
}

\value{

\item{V}{Matrix of dyadic dominance relationships.}
}

\references{

Landau, H. G. (1951). On dominance relations and the structure of animal societies: I. Effect of inherent characteristics. \emph{Bulletin of Mathematical Biophysics}, \emph{13}, 1-19.

Schmid, V. S., & de Vries, H. (2013). Finding a dominance order most consistent with a linear hierarchy: an improved algorithm for the I&SI method. \emph{Animal Behaviour}, \emph{86}, 1097-1105.

de Vries, H. (1995). An improved test of linearity in dominance hierarchies containing unknown or tied relationships. \emph{Animal Behaviour}, \emph{50}, 1375-1389.
}

\author{
David Leiva \email{dleivaur@ub.edu} & Han de Vries \email{J.deVries1@uu.nl}.
}


\examples{

  set.seed(123)
  X <- matrix(floor(runif(20*20,min=0,max=20)),nrow=20,byrow=TRUE)
  diag(X)<-0
  freq2abil(X)
}

\keyword{misc}


