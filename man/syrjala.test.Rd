\encoding{latin1}
\name{syrjala.test}
\alias{syrjala.test}
\alias{plot.ecespa.syrjala}
\alias{print.ecespa.syrjala}

\title{ Syrjala's test for the difference between the spatial distributions of two populations }
\description{
  Computes a two-sample Cramer-von Mises type test for a difference between the spatial distribution of two populations.
  It is designed to be sensitive to differences in the way the populations are distributed across the study area but insensitive to 
  a difference in abundance between the two populations.
}
\usage{
syrjala.test(ppp1, ppp2, nsim = 999)

#S3 method for class 'ecespa.syrjala'
plot.ecespa.syrjala(x, ...)
}

\arguments{
  \item{ppp1}{ A marked point pattern, with the \code{\link[spatstat]{ppp}} format of spatstat, 
  representing the values of some parameter measured on the corresponding sampling locations. }
  \item{ppp2}{ A marked point pattern, with the \code{\link[spatstat]{ppp}} format of spatstat, 
  representing the values of some other parameter measured on the same locations than \code{ppp1}. }
  \item{nsim}{ Number of permutations. }

\item{x}{ An object of class '\code{ecespa.syrjala}' resulting from \code{syrjala.test} }
  \item{...}{ Graphical parameters passed to \code{\link{hist}}.}
}
\details{
  The null hypothesis of  Syrjala's test is that across the study area, 
  the normalized distributions of the two populations are the same (Syrjala, 1996). 
  Population density data are collected at \eqn{K} sampling locations on two populations. 
  Let \eqn{(xk, yk)} denote the coordinates of the \eqn{kth} sampling location (\eqn{ k= 1,...,K} ); 
  let \eqn{d.i(xk, yk)} denote the sample density at the \eqn{Kth} sampling location of the \eqn{ith} population.
   To construct a test that is independent of the population sizes, the observed density data is first normalized:
   \deqn{gamma.i(xk, yk) = di(xk, yk)/Di,} where \eqn{Di} is the sum  of \eqn{d.i(xk, yk)} observations  across the \eqn{K} sampling locations. 
   The value of the cumulative distribution function at the location \eqn{(xk, yk)}  for the \eqn{ith} population, denoted \eqn{GAMMA.i(xk, yk)},
   is the sum of all normalized density observations, \eqn{gamma.i(xk, yk)}, whose location \eqn{(x, y)} is such that 
   \eqn{x <= xk} and \eqn{y <= yk}. The statistic proposed by Syrjala to test the null hypothesis is the square of the difference 
   between the cumulative distribution functions  of the two populations,  summed over all sampling locations,
   that is \deqn{psi = sum{GAMMA.1(xk, yk) - GAMMA.2(xk, yk)}^2.} 
   As \eqn{psi} is not invariant with respect to the 'corner'  of the rectangle enclosing the study area that is chosen as the origin of the coordinate sytem,
   \eqn{psi}  is computed  four times, one with each corner as the origin, and the average \eqn{psi} is employed as the test statistic.
   The level of significance of the observed \eqn{psi} is determined from its position  in the ordered set of test statistic values from all \eqn{2^K} pairwise 
   permutations (that is approximated from a large number of randomly selected permutations).
   
}
\value{
  \code{syrjala.test} returns an object of class '\code{ecespa.syrjala}', basically a list with the following elements: 
  \item{psi.obs }{The observed (averaged) \eqn{psi} statistic.}
  \item{psi.sim }{A vector with the \code{nsim} simulated \eqn{psi}'s statistics.}
  \item{datanames }{A vector with the names of the two point patterns whose spatial congruence is been analyzed.}
  \item{nsim }{The number of permutations employed in the test.}
  The S3 plot method plots an histogram with the distribution of the simulated \eqn{psi}'s statistics and draws the observed \eqn{psi} as a vertical line.
}
\references{


\enc{Rey-Benayas, J.M., de la Montaña, E., Pérez-Camacho, L., de la Cruz, M., Moreno, D., Parejo, J.L. and Suárez-Seoane, S.}{Rey-Benayas, J.M., de la Montaña, E., Perez-Camacho, L., de la Cruz, M., Moreno, D., Parejo, J.L. and Suarez-Seoane, S.} 2008. 
Inter-annual dynamics and spatial congruence of a nocturnal bird assemblage inhabiting a Mediterranean agricultural mosaic. 
\emph{Submitted}.


Syrjala, S. E. 1996. A statistical test for a difference between the spatial distributions of two populations. \emph{Ecology} 77: 75-80. 

}

\author{ Marcelino de la Cruz Rot \email{marcelino.delacruz@upm.es}}
\note{ One can use function \code{\link{haz.ppp}} to easily build the \code{ppp} objects from a \code{data.frame} with only three columns (x-coordinate,
y-coordinate and abundance).  

This function has been employed to compute Syrjala's test in Rey-Benayas et al. (2008).
}
\section{Warning }{The test requires both populations being sampled in exactly the same sampling locations. 
Althoug this implementation employs \code{\link[spatstat]{ppp}}'s as the supporting data format, this kind of data are \bold{not} spatial point patterns.
They cannot be analysed with the usual tools employed for marked point patterns.} 


\examples{
    \dontrun{
   
   data(syr1); data(syr2); data(syr3)
   
   plot(syrjala.test(syr1, syr2, nsim=999)) #non-significant difference of distributions
   
   plot(syrjala.test(syr1, syr3, nsim=999)) #significant difference of distributions
  
  }
}
\keyword{ spatial }
\keyword{ distribution }
