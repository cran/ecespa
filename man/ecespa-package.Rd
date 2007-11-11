\encoding{latin1}
\name{ecespa}
\alias{ecespa}
\alias{ecespa}
\docType{package}
\title{
Functions for spatial point pattern analysis in ecology

}
\description{
Some wrappers, functions and data sets for spatial point pattern analysis, with an ecological bias. 
}
\details{
\tabular{ll}{
Package: \tab ecespa\cr
Type: \tab Package\cr
Version: \tab 1.0-3\cr
Date: \tab 2007-11-11\cr
License: \tab GPL (version 2 or later)\cr
}
}
\author{
Marcelino de la Cruz Rot, with contributions of Philip M. Dixon and heavily borrowing  Baddeley's & Turner's \pkg{spatstat} code.

Mantainer: Marcelino de la Cruz Rot \email{marcelino.delacruz@upm.es}
}

\references{
De la Cruz, M. 2007. \enc{Métodos para analizar datos puntuales. 
  En: \emph{Introducción al Análisis Espacial de Datos en Ecología y Ciencias Ambientales: Métodos y Aplicaciones } 
  (eds. Maestre, F. T., Escudero, A. y Bonet, A.), pp 000-000. Asociación Española de Ecología Terrestre, Universidad Rey Juan Carlos
  y Caja de Ahorros del Mediterráneo, Madrid.}{Metodos para analizar datos puntuales. 
  En: \emph{Introduccion al Analisis Espacial de Datos en Ecologia y Ciencias Ambientales: Metodos y Aplicaciones } 
  (eds. Maestre, F. T., Escudero, A. y Bonet, A.), pp 000-000. Asociacion Española de Ecologia Terrestre, Universidad Rey Juan Carlos
  y Caja de Ahorros del Mediterraneo, Madrid.}
  
De la Cruz, M., Romao, R.L. & Escudero, A. 2007. Where do seedlings go? A spatio-temporal analysis of
 early mortality in a semiarid specialist. \emph{Submitted}.

Diggle, P. J. 2003. \emph{Statistical analysis of spatial point patterns}. Arnold, London. 

Dixon, P.M. 2002. Nearest-neighbor contingency table analysis of spatial segregation for several species.
  \emph{Ecoscience}, \bold{9} (2): 142-151.  
  
Penttinen, A. 2006. Statistics for Marked Point Patterns. In \emph{The Yearbook of the Finnish Statistical Society}, pp. 70-91. 

}


\keyword{ package }
\keyword{ spatial }

\examples{
\dontrun{


###  Summarize the joint pattern of points and marks at different scales
###  with the normalized mark-weighted K-function (Penttinen, 2006). 
###  Compare this function in two consecutive cohorts of Helianthemum
###  squamatum seedlings:

data(seedlings1)

data(seedlings2)
  
s1km <- Kmm(seedlings1, r=1:100)
  
s2km <- Kmm(seedlings2, r=1:100)
  
plot(s1km$r, s1km$Kmm.n, type="l", lty=1, lwd=3, ylim=c(0.6, 1.2),
         xlab="r (cm)", ylab= expression (K[mm](r)),
	 main="Mark-weighted K-function of Hs seedling cohorts")
  
lines(s2km$r, s2km$Kmm.n, lty=2,lwd=3)
  
abline(h=1, lwd=2, lty=3)
  
legend(x=60, y=1.2, legend=c("HsC1", "HsC2", "Ho:"),
	 lty=c(1, 2, 3), lwd=c(3, 2, 2), bty="n")



###  Explore the local relationships between marks and locations (e.g. size 
###  of one cohort of H. squamatum seedlings). Map the marked point pattern 
###  to a random field for visual inspection, with the normalized mark-sum
###  measure (Penttinen, 2006).

data(seedlings1)

seed.m <- marksum(seedlings1, R=20)

marksum.plot(seed.m, what="marksum")  # raw mark-sum measure

marksum.plot(seed.m, what="pointsum") # point sum measure
   
marksum.plot(seed.m,  what="normalized") # normalized  mark-sum measure





###  Test asociation/repulsion between  a "fixed" pattern (e.g. adult
###  H. squamatum plants) and a "variable" pattern (e.g. of surviving and 
###  dead seedlings), with 2.5\% and 97.5\% envelopes of 999 random 
###  labellings (De la Cruz & al. 2007).

data(Helianthemum)


cosa <- K012(Helianthemum, fijo="adultHS", i="deadpl", j="survpl",
             r=seq(0,200,le=201), nsim=999, nrank=25, correction="isotropic")

plot(cosa$k01, sqrt(./pi)-r~r,  col=c(3, 1, 3), lty=c(3, 1, 3), las=1,
         ylab=expression(L[12]), xlim=c(0, 200), 
	 main="adult HS vs. dead seedlings")

plot(cosa$k02, sqrt(./pi)-r~r, col=c(3, 1, 3), lty=c(3, 1, 3), las=1, 
         ylab=expression(L[12]), xlim=c(0, 200),
 	 main="adult HS vs. surviving seedlings")




###  Test differences of agregation and segregation between two patterns, 
###  e.g. surviving and dying H. squamatum seedlings (De la Cruz & al. 2007). 

data(Helianthemum)

cosa12 <- K1K2(Helianthemum, j="deadpl", i="survpl", r=seq(0,200,le=201),
		 nsim=99, nrank=1, correction="isotropic")

plot(cosa12$k1k2, lty=c(2, 1, 2), col=c(2, 1, 2), xlim=c(0, 200),
         main= "survival- death")

plot(cosa12$k1k12, lty=c(2, 1, 2), col=c(2, 1, 2), xlim=c(0, 200),
	 main="segregation of surviving seedlings")

plot(cosa12$k2k12, lty=c(2, 1, 2), col=c(2, 1, 2), xlim=c(0, 200),
         main= "segregation of dying seedlings")




###  Test segregation based on the counts in the contingency table
###  of nearest neighbors in a multitype point pattern (Dixon, 2002)

data(swamp)

dixon2002(swamp,nsim=99)




###  Fit the Poisson cluster point process to a point pattern with the method of
###  minimum contrast (Diggle 2003).

data(gypsophylous)

# Estimate K function ("Kobs").

gyps.env <- envelope(gypsophylous, Kest, correction="iso")

plot(gyps.env, sqrt(./pi)-r~r)

# Fit Poisson Cluster Process. The limits of integration 
# rmin and rmax are setup to 0 and 60, respectively. 

cosa.pc <- pc.estK(Kobs = gyps.env$obs[gyps.env$r<=60],
		           r = gyps.env$r[gyps.env$r<=60])

# Add fitted Kclust function to the plot.

lines(gyps.env$r,sqrt(Kclust(gyps.env$r, cosa.pc$sigma2,cosa.pc$rho)/pi)-gyps.env$r,
       lty=2, lwd=3, col="purple")

# A kind of pointwise test of the gypsophylous pattern been a realisation
# of the fitted model, simulating with sim.poissonc and using function J (Jest).

gyps.env.sim <- envelope(gypsophylous, Jest, 
                    simulate=expression(sim.poissonc(gypsophylous,
		    sigma=sqrt(cosa.pc$sigma2), rho=cosa.pc$rho)))

plot(gyps.env.sim,  main="")


}
}


