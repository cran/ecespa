`syrjala.test` <-
function(ppp1, ppp2, nsim=999){
   
   require(spatstat)
   datanames <- c(deparse(substitute(ppp1)), deparse(substitute(ppp2)))   
   
   # function gamma: computes cummulative distribution function at the location(xk, yk)
   gammaf <- function(cosa.ppp){
              # First, normalize the observed density data
#cosa.ppp$marks <- cosa.ppp$marks/sum(cosa.ppp$marks)
   
      # Second, compute cumulative distribution at each location "k"
gamma <- NULL
for(k in 1: cosa.ppp$n){
        gamma <- c(gamma,
                        sum(cosa.ppp$marks[cosa.ppp$x<=cosa.ppp$x[k] & cosa.ppp$y<=cosa.ppp$y[k]]))
                }
                return(gamma)
        }

    # function psi: computes the squared difference between the two cumulative distribution functions 
    psi <- function(ppp1, ppp2){
gamma1 <- gammaf (ppp1)
gamma2 <- gammaf (ppp2)
psi <- sum( (gamma1-gamma2)^2)
return(psi)
}


    # function psimean: computes the mean psi (i.e among different psi's computed on cummulative distributions
    # relative to the different four corners of the study region). 
    psimean <- function(ppp1, ppp2){

psi1 <- psi(ppp1,ppp2)
psi2 <- psi(affine(ppp1, mat=diag(c(-1,-1))), affine(ppp2, mat=diag(c(-1,-1))))
psi3 <- psi(affine(ppp1, mat=diag(c(1,-1))), affine(ppp2, mat=diag(c(1,-1))))
psi4 <- psi(affine(ppp1, mat=diag(c(-1,1))), affine(ppp2, mat=diag(c(-1,1))))
psimean <- (psi1 + psi2 + psi3 + psi4)/4
return(psimean)
}

      #function rpsi: pairwise permutation of the normalized density observation
      ## OJO: en la implementaci<f3>n original para las nocturnas, lo que se permut<f3> fu<e9> la observaci<f3>n bruta!!      
      # pepesim es una lista de listas (cada elemento lleva los dos point patterns)
      rpsi <- function(pepesim){     
                         marcas <- cbind(ppp1$marks,ppp2$marks)
                         marcas <- t(apply(marcas,1,sample))
                         ppp1$marks <- marcas[,1]
                         ppp2$marks <- marcas[,2]
                       #c<e1>lculo del psi simulado
                       psi.sim <-  psimean (ppp1, ppp2)
        }


# Compute syrjala test and simulations:

# First, normalize the observed density data
ppp1$marks <- ppp1$marks/sum(ppp1$marks)
           ppp2$marks <- ppp2$marks/sum(ppp2$marks)


#Second, compute observed psi
               psi.obs <- psimean(ppp1,ppp2)
   

#Third, instead of the classical loop, generate a list of pairs of populations, 
# permute  pairwisely the density and compute psi sim

     pepesim <- list(pepesim=list(ppp1=ppp1, ppp2=ppp2))
     pepesim <- rep(pepesim, nsim) 
     psi.sim <- sapply(pepesim, rpsi)
     
       

   result <- list(psi.obs=psi.obs, psi.sim=psi.sim, datanames=datanames, nsim=nsim)
   class(result) <- c("ecespa.syrjala", class(result))
   return(result)


}

