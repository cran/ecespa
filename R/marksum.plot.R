`marksum.plot` <-
function(mimarksum, what="normalized", gris=FALSE,
                      contour=FALSE, interpol=TRUE, leip=100,
                      main="",...){
   
   require (akima)
   whatn=NULL
   if(what=="normalized") whatn=1
   if(what=="marksum") whatn=2
   if(what=="pointsum") whatn=3
   colorido=heat.colors(12)
   if(gris==TRUE) colorido=gray((0:64)/64)
   
   if (interpol==FALSE){
       image(sort(unique(mimarksum$grid$x)),sort(unique(mimarksum$grid$y)),
             matrix(mimarksum[[whatn]], mimarksum$nx, mimarksum$ny),
             xlab="",ylab="",col=colorido, main=main)
       if (contour==TRUE){
       contour(sort(unique(mimarksum$grid$x)),sort(unique(mimarksum$grid$y)),
               matrix(mimarksum[[whatn]], mimarksum$nx, mimarksum$ny),
               add=TRUE, ...)
       }

   }
   if (interpol==TRUE){
       cosainterp = interp(mimarksum$grid$x, mimarksum$grid$y,
                           mimarksum[[whatn]], 
                           xo=seq(mimarksum$window$x[1], mimarksum$window$x[2],le=leip),
                           yo=seq(mimarksum$window$y[1], mimarksum$window$y[2],le=leip),
                            )   
        image(cosainterp, xlab="", ylab="", col=colorido,main=main)
        if (contour==TRUE){
            contour(cosainterp,add=TRUE,...)
        }
   }
}

