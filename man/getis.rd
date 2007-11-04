\name{getis}
\alias{getis}
\alias{getis.plot}
\title{ Neighbourhood density function }
\description{
  Computes and plots the neighbourhood density function, a local version of the K-function defined by Getis and Franklin (1987). 
}
\usage{
getis(mippp, nx = 30, ny = 30, R = 10)
getis.plot(getis.obj, type="k", interp=100, color=tim.colors(64), contour=TRUE, points=TRUE)
}
\arguments{
  \item{mippp}{A point pattern. An object with the \code{\link[spatstat]{ppp}} format of \code{spatstat}. }
  \item{nx}{Grid dimensions (for estimation) in the x-side. }
  \item{ny}{Grid dimensions (for estimation) in the y-side. }
  \item{R}{Radius. The distance argument \emph{ r} at which the function K should be computed. }
  
  \item{getis.obj}{Result of applying \code{getis} to a point pattern. }
  \item{type}{Type of local statistics to be ploted. One of \code{k} (local-K), \code{l} (local-L), \code{n} (local-n) or \code{d} (deviations from CSR).}
  \item{interp}{Number of points in the side of the grid of points to interpolate the results.}
  \item{color}{A list of colors such as that generated by \code{\link{rainbow}}, \code{\link{heat.colors}}, \code{\link{topo.colors}}, \code{\link{terrain.colors}} or similar functions.}
  \item{contour}{logical; if TRUE, add a contour to current plot}
  \item{points}{logical; if TRUE, add the point pattern to current plot}
}
\details{
  Getis and Franklin (1987) proposed the neigbourhood density function, a local version of Ripley's L- function.
 Given a spatial point pattern X, the neigbourhood density function associated with the \emph{i}th point in \code{X} is computed by

\eqn{L[i](r) = sqrt((a/((n-1))*pi))*sum[j]e[i,j])}
where the sum is over all points \emph{ j != i} that lie within a distance \emph{r} of the \emph{i}th point, \emph{a} is the area of the observation window,
\emph{n} is the number of points in \code{X}, and \emph{e[i,j]} is the isotropic edge correction term (as described in \code{\link[spatstat]{Kest}}). The value of \emph{L[i](r)} can also
 be interpreted as one of the summands that contributes to the global estimate of the L function. 

The command \code{getis} actually computes the local K-function using \code{\link[spatstat]{Kcross}}. As the main objective of \code{getis} is to map the local density function,  
as sugested by Gestis and Franklin (1987: 476) a grid of  points (whose density is controled by \code{nx} and \code{ny}),  is used to accurately estimate the
functions in empty or sparse areas. The command \code{\link{getis.plot}}  plots the spatial distribution of  the local K function or other related local statistics, such as 
\eqn{L[i](r)},   \eqn{n[i](r)} [=\eqn{ lambda*K[i](r)}] or the deviations from  the expected value of  local  L  under CSR [= \eqn{L[i](r) -r}].  It uses the function 
\code{\link[akima]{interp}} in \code{akima} package to interpolate the results.

}
\value{
  \code{getis} gives a list with the following elements:
  \item{x }{x coordinates of pattern points (ahead) and grid points}
  \item{y }{y coordinates of pattern points (ahead) and grid points}
  \item{klocal }{estimate of local \eqn{K[i](r)} at the point pattern points}
  \item{klocalgrid }{estimate of local \eqn{K[i](r)} at the grid points}
  \item{R }{distance \eqn{r} at which the estimation is made}
  \item{ppp }{original point pattern}

  \code{getis.plot} plots an interpolated map of the selected local statistics
}
\note{
As \code{getis.plot} interpolates over rectangular grid of points, it is not apropriate to map irregular windows. In those cases, \code{\link[spatstat]{smooth.ppp}} of \code{spatstat}
can be used to interpolate the local statistics (see examples)
}
\references{ Getis, A. and Franklin, J. (1987) Second-order neighbourhood analysis of mapped point patterns. \emph{Ecology} \bold{68}: 473-477
 }
\author{ Marcelino de la Cruz Rot \email{marcelino.delacruz@upm.es} }
\seealso{ \code{\link[spatstat]{localK}}, a different approach in \code{spatstat} }
\examples{
 \dontrun{
  ## Compare with fig. 5b of Getis & Franklin (1987: 476):
  data(ponderosa)
  ponderosa12=getis(ponderosa, nx = 30, ny = 30, R = 12)
  getis.plot(ponderosa12, type = "l")

  ## Plot the same, using smooth.ppp in spatstat
  ponderosa.12=setmarks(ponderosa, ponderosa12$klocal)
  Z <- smooth.ppp(ponderosa.12, sigma=5, dimyx=256)
  plot(Z, col=topo.colors(128), main="smoothed neighbourhood density")
  contour(Z, add=TRUE)
  points(ponderosa, pch=16, cex=0.5) 
  
  ## Example with irregular window:
  data(letterR)
  X <- rpoispp(50, win=letterR)
  X.g <- getis(X, R=0.2)
  X2 <- setmarks(X, X.g$klocal)
  Z <- smooth.ppp(X2, sigma=0.05, dimxy=256)
  
  plot(Z, col=topo.colors(128), main="smoothed neighbourhood density")
  contour(Z, add=TRUE)
  points(X, pch=16, cex=0.5)
  
  
  
    }
}
\keyword{ spatial }