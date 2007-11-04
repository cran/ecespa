`marksum` <-
function(mippp, R=10, nx=30, ny=30){
   require(spatstat)
   verifyclass(mippp, "ppp")
   if(mippp$window$type!="rectangle") stop ("marksum only implemented for rectangular windows")  
   if(is.marked(mippp)!=T) stop ("marksum only implemented for **marked** patterns")
   
  ##meter stopif para ppp, marked, window=rectangle, etc

   grid=gridcenters(mippp$window,nx,ny)

   grid.ppp= ppp(x=grid$x,y=grid$y,marks=rep(0,length(grid$x)),window=mippp$window)
   
   prueba = superimpose(grid.ppp, mippp)
   marksum = markstat(prueba, sum, R=R)
   marksum = marksum[1:length(grid$x)] ## we are only interested in the grid stats
   
   ##now we count all the points within a distance R of the grid points
   pointsum = markstat(prueba, length, R=R) # this will count both grid and mippp points
   pointsum = pointsum[1:length(grid$x)] ## we are only interested in the grid stats
   
## correction of excess points in pointsum. We must substract the grid points summed
## in the previous markstat
  minus = markstat(grid.ppp, length, R=R) # this will count exclusively the grid points
  pointsum = pointsum-minus
  normalized.marksum=marksum/pointsum  
  normalized.marksum[marksum==0]=0 ## correction div by 0
  return(list(normalized=normalized.marksum, marksum=marksum,pointsum=pointsum,
              minus = minus,grid=grid,nx=nx,ny=ny,R=R,window=mippp$window))

}

