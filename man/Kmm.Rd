\name{Kmm}
\alias{Kmm}
\title{ Mark-weighted K-function }
\description{
  This is a functional data summary for marked point patterns that measures de joint pattern of points and marks at different scales determined by \eqn{r}.
}
\usage{
Kmm(mippp, r = 1:100)
}
\arguments{
   \item{mippp}{ A marked point pattern. An object with the \code{\link[spatstat]{ppp}} format of \pkg{spatstat}. }
  \item{r}{ Sequence of distances at which \code{Kmm} is estimated }
}
\details{
  The mark-weighted K-function  of a stationary marked point process \eqn{X} 
  is defined so that \deqn{lambda*Kmm(r) = Eo[sum(mo*mn)]/mu^2} where \eqn{lambda} is the intensity of the process, 
  i.e. the expected number of points of X per unit area, \eqn{Eo[ ] }denotes expectation (given that there is a point at the origin); 
  \eqn{m0} and \eqn{mn} are the marks attached to every two points of the process separated by a distance \eqn{<= r}  and \eqn{mu}
  is the mean mark. It measures the joint pattern of marks and points at the scales determmined by \eqn{r}. If all the marks are set to 1,
  then \eqn{lambda*Kmm(r)} equals the expected number of additional random points within 
  a distance \eqn{r} of a typical random point of \eqn{X}, i.e. \eqn{Kmm} becomes the conventional Ripley's \eqn{K}-function for unmarked point processes.
  As the \eqn{K}-function measures clustering or regularity among the points regardless of the marks, one can separate clustering of marks with the
  \emph{normalized weighted K-function} \deqn{bar(Kmm) = Kmm(r)/K(r)} If the process is independently marked, \eqn{Kmm(r)} equals \eqn{K(r)} so the 
  normalized mark-weighted K-function will equal 1 for all distances \eqn{r}.
  }
\value{
   \item{r }{Sequence of distances at which \code{Kmm} is estimated }
  \item{Kmm }{Mark-weighted K-function}
  \item{Kmm.n }{Normalized mark-weighted K-function}
}
\references{ 
Penttinen, A. 2006. Statistics for Marked Point Patterns. In \emph{The Yearbook of the Finnish Statistical Society}, pp. 70-91. 
}
\author{ Marcelino de la Cruz Rot \email{marcelino.delacruz@upm.es}}

\note{This implementation estimates \eqn{Kmm(r)} without any correction of border effects, so it must be used with caution. However, as \eqn{K(r)} is 
also estimed without correction it migth compensate the border effects on the normalized Kmm-function.
}
\seealso{\code{\link[spatstat]{markcorr}} }
\examples{
  \dontrun{
  data(seedlings1)
  data(seedlings2)
  s1km <- Kmm(seedlings1, r=1:100)
  s2km <- Kmm(seedlings2, r=1:100)
  plot(s1km$r, s1km$Kmm.n, type="l", lty=1, lwd=3, ylim=c(0.6, 1.2),
         xlab="r (cm)", ylab= expression (K[mm](r)),
	 main="Mark-weighted K-function of Hs seedling cohorts")
  lines(s2km$r, s2km$Kmm.n, lty=2,lwd=3)
  abline(h=1, lwd=2, lty=3)
  legend(x=60, y=1.2, legend=c("HsC1","HsC2","Ho:"),
	 lty=c(1,2,3), lwd=c(3,2,2), bty="n")
  }
}
\keyword{spatial }

