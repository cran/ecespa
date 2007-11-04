`Kmm` <-
function(mippp, r=1:100){
   mippp1 = setmarks(mippp,1)
   lambda = mippp$n/area.owin(mippp$window)
   mu=mean(mippp$marks)
   Kmm=NULL
   Kmm1=NULL
   m0 = markstat(mippp,sum,R=0)
   m01 = markstat(mippp1,sum,R=0)

   for ( i in 1:length(r)){
     progressreport(i, length(r))
     mr = markstat(mippp,sum,R=r[i])
     mr1 = markstat(mippp1,sum,R=r[i])
     sumatorio=(mr-m0)*m0
     sumatorio1=(mr1-m01)*m01
     E0 = mean(sumatorio)
     E01 = mean(sumatorio1)
     Kmm = c(Kmm,E0/(lambda*mu^2))
     Kmm1 = c(Kmm1,E01/(lambda)) ## Función K de puntos con marca = 1
   }
   return(list(r=r, Kmm=Kmm, Kmm.n=Kmm/Kmm1))
}

