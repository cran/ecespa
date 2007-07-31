`plot.ecespa.marksum` <-
function(x, what="normalized", contour=F, grid=F, ...)
{
 require (spatstat)
 if(what=="normalized") {cosa <- x$normalized; what="normalized mark-sum"}
 if(what=="pointsum") {cosa <-  x$pointsum; what="point-sum"}
 if(what=="marksum") {cosa <-  x$marksum; what="mark-sum"}
 
 plot(smooth.ppp(setmarks(x$grid.ppp, cosa),...), main="")
 title(main=paste(x$dataname,"\n",noquote(what), "measure; R=", x$R))
 if (contour==TRUE) contour(smooth.ppp(setmarks(x$grid.ppp, cosa),...), add=T)
 if (grid==TRUE) plot(setmarks(x$grid.ppp, cosa), add=T)
}

