\name{gypsophylous}
\alias{gypsophylous}
\docType{data}
\title{ Spatial point pattern of a plant community}
\description{
  Locations of  plants  in a gypsophylous plant community in Central Spain. These are part of the data collected by Romao (2003) that have been
 analyzed several times (Escudero \emph{et al.} 2005, De la Cruz 2006).
 The coordinates of the plans are given in cm. 
}
\usage{data(gypsophylous)}
\format{
 An object of class "ppp"  of \code{spatstat} representing the point pattern of plants locations. See \code{\link[spatstat]{ppp.object}} 
  for details of the format.
}
\source{
Romao, R.L. 2003. \emph{Estructura espacial de comunidades de gipsofitos: interacciones bioticas
 y constricciones abioticas.} Tesis Doctoral Inedita. Universidad Politecnica de Madrid.
}
\references{
  De la Cruz, M. 2006. Introduccion al analisis de datos mapeados o algunas de las (muchas) cosas
 que puedo hacer si tengo coordenadas. \emph{Ecosistemas}. 2006/3. 
\url{http://www.revistaecosistemas.net/articulo.asp?Id=488&Id_Categoria=1&tipo=portada}.

Escudero, A., Romao, R.L., De la Cruz, M. & Maestre, F. 2005. Spatial pattern and neighbour effects on 
\emph{Helianthemum squamatum} seedlings in a Mediterranean gypsum community. \emph{J. Veg. Sci.}, \bold{16}: 383-390.

}
\examples{
data(gypsophylous)
plot(gypsophylous)
}
\keyword{datasets}
