`Kci` <-
function(mod1,mod2,correction="trans", nsim=99,
                         ngrid=200, nrep=1e5, r=NULL, simu="both", 
                         spctype=1)
{
   verbose <- FALSE # 
   ## datos b<e1>sicos
   modnamea <- deparse(substitute(mod1))
   modnameb <- deparse(substitute(mod2))
   
   if(inherits(mod1,"ppm")){  
      I.ppp <- mod1$Q$data
      lambdaI <- predict(ppm(mod1$Q$data, mod1$trend), type="trend", ngrid=ngrid)
      Isim <- "ppm"
      dataname.a <- mod1$call[[2]]
   }
   if(inherits(mod1,"ecespa.minconfit")){  
      I.ppp <- mod1$data
      lambdaI <- mod1$lambda
      Isim <- "spc"
      dataname.a <- mod1$dataname
   }

   if(inherits(mod2,"ppm")){  
      J.ppp <- mod2$Q$data
      lambdaJ <- predict(ppm(mod2$Q$data, mod2$trend), type="trend", ngrid=ngrid)
      Jsim <- "ppm"
      dataname.b <- mod2$call[[2]]
   }
   if(inherits(mod2,"ecespa.minconfit")){  
      J.ppp <- mod2$data
      lambdaJ <- mod2$lambda
      Jsim <- "spc"
      dataname.b <- mod2$dataname
   }

   ## C<e1>lculo de los Ki de cada patr<f3>n
   Kia <- Kinhom(I.ppp, lambdaI, correction=correction, r=r)
   mi.r <- Kia$r
   Kia <- Kia[[3]]
   Kib <- Kinhom(J.ppp, lambdaJ, correction=correction, r = mi.r)
   Kib <- Kib[[3]]
   
   ## generaci<f3>n del patr<f3>n multivariado y c<e1>lculo de Kci observada
   ## ojo: por el borde y la inhomogeneidad no es lo mismo Kab que Kba
   
   IJ.ppp <- superimpose(a=I.ppp, b=J.ppp)   
   Kci.ab.o <- Kcross.inhom(IJ.ppp, i="a", j="b", lambdaI, lambdaJ,
                            correction= correction, r=mi.r)[[3]]
   Kci.ba.o <- Kcross.inhom(IJ.ppp, i="b", j="a", lambdaJ, lambdaI,
                            correction= correction, r=mi.r)[[3]]

   Kia.s <- NULL ##para la K univariada del patr<f3>n que se simula
   Kib.s <- NULL
   Kci.ab.s <- NULL
   Kci.ba.s <- NULL
   # Isim.ppp <- I.ppp
   ## start simulations
   for (i in 1: nsim){ 
      progressreport(i,nsim)
      ## simulate from second model:
      if (Jsim=="ppm"){     Jsim.ppp <- rmh(mod2, start=list(x.start=J.ppp),
                                            control=list(p=1, nrep=nrep),
                                             verbose=verbose)}
      else if (Jsim=="spc") Jsim.ppp <- rIPCP (mod2, type=spctype)
      
      ## aseguramos que no haya NAs en el vector simulado de lambdas
      dentro <- !is.na(lambdaJ[Jsim.ppp, drop=FALSE]) 
      Jsim.ppp <- Jsim.ppp[dentro] 

      ## simulated multivariate PP
      IJs.ppp <- superimpose(a=I.ppp, b=Jsim.ppp, W=I.ppp$w)
      IsJ.ppp <- IJs.ppp # si solo se simula J, creamos dos multipatrones iguales

      if(simu=="both"){
          if(Isim=="ppm"){    Isim.ppp <- rmh(mod1, start=list(x.start=I.ppp),
                                               control=list(p=1, nrep=nrep),
                                               verbose=verbose)}
          else if (Isim=="spc") Isim.ppp <- rIPCP (mod1, type=spctype)
          
          ## aseguramos que no haya NAs en el vector simulado de lambdas
          dentro <- !is.na(lambdaI[Isim.ppp, drop=FALSE]) 
          Isim.ppp <- Isim.ppp[dentro]

          ## simulated multivariate PP:
          IsJ.ppp <- superimpose(a=Isim.ppp, b=J.ppp, W=I.ppp$w)
      } 

      ## K simuladas
      Kib.s <- cbind(Kib.s, Kinhom(Jsim.ppp, lambdaJ, correction=correction,
                                  r = mi.r, nlarge=Inf)[[3]])
      Kci.ab.s <- cbind(Kci.ab.s, Kcross.inhom(IJs.ppp, i="a", j="b",
                                     lambdaI, lambdaJ, r=mi.r,
                                     correction=correction)[[3]])
       Kci.ba.s <- cbind(Kci.ba.s, Kcross.inhom(IsJ.ppp, i="b", j="a",
                                     lambdaJ, lambdaI, r=mi.r,
                                     correction=correction)[[3]])
       if(simu=="both"){Kia.s <- cbind(Kia.s, Kinhom(Isim.ppp, lambdaI, 
                                                 correction=correction, r = mi.r, 
                                                 nlarge=Inf)[[3]])}
   }
   

   result <- list(r=mi.r, kia = Kia, kib=Kib, kci.ab.o=Kci.ab.o,
                kci.ba.o=Kci.ba.o, kci.ab.s=Kci.ab.s, kci.ba.s=Kci.ba.s,
                kib.s=Kib.s, kia.s=Kia.s, datanamea=dataname.a, datanameb=dataname.b,
modnamea=modnamea, modnameb=modnameb, type="Kci")
   class(result) <- c("ecespa.kci", class(result))
   return(result)
}

