\name{Helianthemum}
\alias{Helianthemum}
\docType{data}
\title{ Spatial point pattern of Helianthemum squamatum adult plants and seedlings }
\description{
  Locations of \emph{H. squamatum} adult plants and seedlings in a 6 m x 7 m plot over gypsum soil 
in Chinchon (near Madrid, Spain). These are part of the data collected by Romao (2003) that have been
 analyzed several times (Escudero \emph{et al.}2005, De la Cruz 2006, De la Cruz ,\emph{et al.} 2007).
 The coordinates of the plans are given in cm. 

}
\usage{data(Helianthemum)}
\format{
 An object of class "ppp"  of \code{spatstat} representing the point pattern of plants locations marked by their type. See \code{\link[spatstat]{ppp.object}} 
  for details of the format.The dataset has 866 points with the following levels: 
  \describe{
    \item{\code{adultHS}}{adult \emph{H. squamatum} plants}
    \item{\code{deadpl}}{dying \emph{H. squamatum} seedlings}
    \item{\code{survpl}}{surviving \emph{H. squamatum} seedlings}
  }
}
\source{
Romao, R.L. 2003. \emph{Estructura espacial de comunidades de gipsofitos: interacciones bioticas
 y constricciones abioticas.} Tesis Doctoral Inedita. Universidad Politecnica de Madrid.
}
\references{
  De la Cruz, M. 2006. Introduccion al analisis de datos mapeados o algunas de las (muchas) cosas
 que puedo hacer si tengo coordenadas. \emph{Ecosistemas}. 2006/3. 
\url{http://www.revistaecosistemas.net/articulo.asp?Id=488&Id_Categoria=1&tipo=portada}.

De la Cruz, M., Romao, R.L. & Escudero, A. 2007. Where do seedlings go? A spatio-temporal analysis of
early mortality in a semiarid specialist. \emph{Ecography, in revision}.

Escudero, A., Romao, R.L., De la Cruz, M. & Maestre, F. 2005. Spatial pattern and neighbour effects on 
\emph{Helianthemum squamatum} seedlings in a Mediterranean gypsum community. \emph{J. Veg. Sci.}, \bold{16}: 383-390.

}
\examples{
data(Helianthemum)
plot(Helianthemum)
}
\keyword{datasets}
