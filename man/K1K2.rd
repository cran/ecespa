\name{K1K2}
\alias{K1K2}
\title{ Differences between univariate and bivariate K-functions }
\description{
  Given two point patterns I and J, \code{K1K2}computes the differences between 
both univariate K-functions (i.e. \eqn{Ki(r)-Kj(r)}) as well as the differences between 
the univariate and the combined bivariate K-function (i.e. \eqn{Ki(r)-K*ij(r)} and \eqn{Kj(r)-K*ij(r)}).
 It also computes simulation envelopes to test  that that the observed differences are within the 
range expected asuming the random labelling hypothesis .
}
\usage{
K1K2(X, i, j, nsim = 99, nrank = 1, r = NULL, correction = "isotropic")
}
\arguments{
  \item{X}{ Multitype marked point pattern. An object with the \code{\link[spatstat]{ppp}} format of \pkg{spatstat}.  }
  \item{i}{ Number or character string identifying the mark value of the  I pattern in X }
  \item{j}{ Number or character string identifying the mark value of the  J pattern in X }
  \item{nsim}{ Number of simulated point patterns to be generated when computing the envelopes.}
  \item{nrank}{ Integer. Rank of the envelope value amongst the nsim simulated values. 
A rank of 1 means that the minimum and maximum simulated values will be used. }
 \item{r}{ Numeric vector. The values of the argument r at which the K functions  should be evaluated. }
  \item{correction}{ A character item selecting any of the options "border", "bord.modif", "isotropic", "Ripley" or
 "translate". It specifies the edge correction(s) to be applied. }
}
\details{
 The indiscriminate use of the raw bivariate functions (mainly the K or the L-bivariate functions) in ecological studies for testing the association/ repulsion
 between different point patterns waste some of the most interesting properties of the K function. One of them is that under the random labelling hypothesis
 every individual pattern would be a random thinning of the corresponding bivariate pattern and therefore \eqn{Ki(r)=Kj(r)= Kij(r)=pi r^2} (Diggle 2003). 
Dixon (2002) sugested that some differences of these functions could provide provide interesting ecological information. For example, \eqn{Ki(r)-Kj(r)},
 has an expected value of 0 for all \eqn{r} distances under random labelling and evaluates the differences in the intensity of aggregation of the two point patterns
 (e.g., in the example bellow, the pattern of drought and herbivory deaths). Other relevant function is \eqn{Ki(r)-Kij(r}) and the complementary \eqn{Kj(r)-Kij(r})
 which evaluate the degree of segregation of every individual pattern, i.e. if every point of the pattern is more -or less- surrounded by other points of the same type
 than would be expected under the random labelling hypothesis. \code{K1K2} uses \eqn{K*ij(r)}, the combined estimator of Upton and Fingleton (a weigthed mean of 
 \eqn{Kij(r)}  and \eqn{Kji(r)}) as computed by \code{\link{Kmulti.uf}}.
}
\value{
  A list with three elements.
  \item{k1k2 }{Difference between \eqn{Ki(r)} and \eqn{Kj(r)}, with simulation envelopes}
  \item{k1k12 }{Difference between \eqn{Ki(r)} and \eqn{K*ij(r)}, with simulation envelopes}
  \item{k2k12 }{Difference between \eqn{Kj(r)} and \eqn{K*ij(r)}, with simulation envelopes}
  \item{}{}
    \item{}{}
\item{}{Each of the above elements is a \code{\link[spatstat]{fv.object}}, essentially a \code{data.frame} with the following items:}
\item{r }{the values of the argument r at which the functions kave been estimated}
\item{hi }{upper envelope of simulations}
\item{lo }{lower envelope of simulations}
together with the observed difference in each case (respectively \code{K1-K2}, \code{K1-K12} and \code{K2-K12}). 
}

\references{ 
 De la Cruz, M. 2006. Introduccion al analisis de datos mapeados o algunas de las (muchas) cosas
 que puedo hacer si tengo coordenadas. \emph{Ecosistemas}. 2006/3. 
\url{http://www.revistaecosistemas.net/articulo.asp?Id=488&Id_Categoria=1&tipo=portada}

De la Cruz, M., Romao, R.L. & Escudero, A. 2007. Where do seedlings go? A spatio-temporal analysis of
early mortality in a semiarid specialist. \emph{Submited}.

Diggle, P.J. 2003. \emph{Statistical analysis of spatial point patterns}. Arnold, London.

Dixon, P. M. 2002. Ripley's K function. In \emph{The encyclopedia of environmetrics} 
(eds. El-Shaarawi, A.H. & Piergorsch, W.W.), pp. 1976-1803. John Wiley & Sons Ltd, NY.

}
\author{ Marcelino de la Cruz \email{marcelino.delacruz@upm.es} }
\examples{
\dontrun{
data(Helianthemum)
cosa12= K1K2(Helianthemum, j="deadpl", i="survpl",r=seq(0,200,le=201),nsim=99,nrank=1,correction="isotropic")

%%## plots of figure 9 in De la Cruz (2006)
plot(cosa12$k1k2,lty=c(2,1,2), col=c(2,1,2),xlim=c(0,200), main="survival- death")

plot(cosa12$k1k12,lty=c(2,1,2), col=c(2,1,2),xlim=c(0,200), main="segregation of surviving seedlings")

plot(cosa12$k2k12,lty=c(2,1,2), col=c(2,1,2),xlim=c(0,200), main="segregation of dying seedlings")}
}

\keyword{ spatial}

