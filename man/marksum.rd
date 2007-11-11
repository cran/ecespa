\name{marksum}
\alias{marksum}
\alias{marksum.plot}
\title{ Mark-sum measure }
\description{
  An exploratory data analysis technique for marked point patterns. The marked point pattern is mapped to a random field for visual inspection.
}
\usage{
marksum(mippp, R = 10, nx = 30, ny = 30)
marksum.plot(mimarksum, what="normalized", gris=FALSE,
            contour=FALSE, interpol=TRUE, leip=100, main="", ...)

}
\arguments{
  \item{mippp}{ A marked point pattern. An object with the \code{\link[spatstat]{ppp}} format of \pkg{spatstat}. }
  \item{R}{ Radius. The distance argument \emph{ r} at which the mark-sum measure should be computed }
  \item{nx}{Grid density (for estimation) in the x-side. }
  \item{ny}{ Grid density (for estimation) in the y-side. }
  
  \item{mimarksum}{Result of applying \code{marksum} to a point pattern.}
  \item{what}{What to plot. One of \code{"marksum"} (raw mark sum measure), \code{"point"} (point sum measure) or \code{"normalized"} (normalized sum measure).} 
  \item{gris}{Logical; if \code{"TRUE"} display map in grey levels.}
  \item{contour}{Logical; if \code{"TRUE"} add contour to map.}
  \item{interpol}{Logical; if \code{"TRUE"} display interpolated results (with function \code{\link[akima]{interp}} of \pkg{akima}).}
  \item{leip}{Number of points in the side of the grid of points to interpolate the results.}
  \item{main}{Text or expression to be displayed as a title in the map.}
  \item{...}{Additional parametrs to \code{\link{contour}}.}
}
\details{
  Penttinen (2006) defines the \emph{mark-sum measure} as a smoothed summary measuring locally the contribution of points and marks. For any fixed location \eqn{x} within the 
  observational   window and a distance \eqn{R}, the mark-sum measure \eqn{ S[R](x)} equals the sum of the marks of the points within the circle of radius
  \eqn{R} with centre in  \eqn{x}. The \emph{point-sum measure} \eqn{ I[R](x)} is defined by him as the sum of points within the circle of radius \eqn{R} with centre
  in  \eqn{x}, and describes the contribution of points locally near \eqn{x}. The \emph{normalized mark-sum measure} describes the contribution of marks 
  near \eqn{x} and is defined (Penttinen, 2006) as
  \deqn{ S.normalized[R](x) = S[R](x)/I[R](x)}
  This implementation of \code{marksum} estimates the mark-sum and the point-sum measures in a grid of points whose density is defined by \code{nx} and 
  \code{ny}.
}
\value{
  \code{marksum} gives a list with the following elements:
   \item{normalized }{Normalized mark-sum measure estimated in the grid points. } 
  \item{marksum }{Raw mark-sum measure estimated in the grid points. } 
  \item{pointsum }{Point-sum measure estimated in the grid points. } 
  \item{minus}{Point-sum of the grid points. For advanced use only.} 
  \item{grid}{ Grid of points. } 
  \item{nx }{Density of the estimating grid  in the x-side. } 
  \item{ny }{Density of the estimating grid  in the x-side. } 
  \item{R}{ Radius. The distance argument \emph{r} at which the mark-sum measure has been computed. }
  \item{window}{Window of the point pattern.}
  
  \code{marksum.plot} plots the selected mark-sum measure.
}
\references{ 
Penttinen, A. 2006. Statistics for Marked Point Patterns. In \emph{The Yearbook of the Finnish Statistical Society}, pp. 70-91. 
}
\author{ Marcelino de la Cruz Rot \email{marcelino.delacruz@upm.es}}

\seealso{ \code{\link{getis}}, related to the point-sum measure, and  \code{\link[spatstat]{markstats}} for designing different implementations }
\examples{
 \dontrun{
   
 data(seedlings1)
   
 seed.m <- marksum(seedlings1, R=20)

 marksum.plot(seed.m, what="marksum")  # raw mark-sum measure

 marksum.plot(seed.m, what="pointsum") # point sum measure
   
 marksum.plot(seed.m,  what="normalized") # normalized  mark-sum measure
  
  }
}  
\keyword{ spatial }
