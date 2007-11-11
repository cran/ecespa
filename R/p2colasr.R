`p2colasr` <-
function(Z,nsim=length(Z)){
   ## calcula el p-value para "two-sided hypotesis test" siguiendo el
   ## m<e9>todo de Agresti & Min (2001) seg<fa>n la descripci<f3>n de Dixon (2002:145)
   ## Z es un vector con el estad<ed>stico observado en la primera posici<f3>n y 
   ## a continuaci<f3>n todos los valores simulados

   if(Z[1]<0){
      p1=rank(Z)[1]/(nsim+1)
      ptodos= 1-(rank(Z)/ (nsim+1))
      if (sum(ptodos[ptodos<p1])==0) p2=p1 else p2=p1+ max(ptodos[ptodos<p1])
   }
   else
    { p1=1-rank(Z)[1]/ (nsim+1)
      ptodos= rank(Z)/(nsim+1)
      if (sum(ptodos[ptodos<p1])==0) p2=p1 else p2=p1+ max(ptodos[ptodos<p1])
   }
return(p2)
}

