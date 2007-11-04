`getis.plot` <-
function(getis.obj, type="k", interp=100, color=tim.colors(64), 
                    contour=TRUE , points=TRUE){
   require(fields)
   require (akima)
   require (spatstat)
   lambda = getis.obj$ppp$n/area.owin(getis.obj$ppp$window)
   if (type=="k") zg = getis.obj$klocalgrid
   if (type== "l") zg = sqrt(getis.obj$klocalgrid/pi)
   if (type== "n") zg = getis.obj$klocalgrid*lambda
   if (type== "d") zg = sqrt(getis.obj$klocalgrid/pi)-getis.obj$R
   seqx= seq(round(getis.obj$ppp$window$xrange)[1],
             round(getis.obj$ppp$window$xrange)[2],
             length=interp)
   seqy= seq(round(getis.obj$ppp$window$yrange)[1],
             round(getis.obj$ppp$window$yrange)[2],
             length=interp)
   image.plot(interp(x=getis.obj$x, y=getis.obj$y, z= zg,
                     xo=seqx, yo=seqy), col=color)
   if(contour==TRUE) contour(interp(x=getis.obj$x, y=getis.obj$y, z= zg,
                                     xo=seqx, yo=seqy), add=TRUE)
   if(points==TRUE) plot(getis.obj$ppp, pch=16,, cex=0.7, add=T)
}

