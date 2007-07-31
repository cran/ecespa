`plot.ecespa.getis` <-
function(x, type="k", interp=100, color=tim.colors(64), 
                    contour=TRUE , points=TRUE,...){
   require(fields)
   require (akima)
   require (spatstat)
   lambda = x$ppp$n/area.owin(x$ppp$window)
   if (type=="k") zg = x$klocalgrid
   if (type== "l") zg = sqrt(x$klocalgrid/pi)
   if (type== "n") zg = x$klocalgrid*lambda
   if (type== "d") zg = sqrt(x$klocalgrid/pi)-x$R
   seqx= seq(round(x$ppp$window$xrange)[1],
             round(x$ppp$window$xrange)[2],
             length=interp)
   seqy= seq(round(x$ppp$window$yrange)[1],
             round(x$ppp$window$yrange)[2],
             length=interp)
   image.plot(interp(x=x$x, y=x$y, z= zg,
                     xo=seqx, yo=seqy), col=color,...)
   if(contour==TRUE) contour(interp(x=x$x, y=x$y, z= zg,
                                     xo=seqx, yo=seqy), add=TRUE)
   if(points==TRUE) plot(x$ppp, pch=16,, cex=0.7, add=T)
}

