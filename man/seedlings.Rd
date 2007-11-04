\name{seedlings}
\alias{seedlings}
\alias{seedlings1}
\alias{seedlings2}
\docType{data}
\title{ Cohorts of Helianthemum squamatum seedlings }
\description{
  Marked point patterns of two consecutive cohorts of seedlings of \emph{H. squamatum} growing in a gypsophylous plant community in Central Spain. The datasets contains the
  locations of the seedlings marked with their heigth. Both the coordinates and the heigth of the seedlings are given in cm.
}
\usage{
data(seedlings1)
data(seedlings2)
}
\format{
  seedlings1 and seedlings2 are objects of class "ppp" representing the point pattern of seedling locations marked by their heights. See \code{\link[spatstat]{ppp.object}} 
  for details of the format. }

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
\dontrun{
data(seedlings1)
 plot(seedlings1)
}
}
\keyword{datasets}
