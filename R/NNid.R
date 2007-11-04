`NNid` <-
function(xy, splancs=TRUE) {
 
# find id of nearest neighbor, given xy positions as a matrix
# uses either of two functions,
# depending on whether splancs or spatialstats is available

# if group is non null, then finds the neighbor that is withing
#   the same group as the point.  Typical use is restricting
#   neighbors to same side of the street.

  if (splancs) { 
     nndistG(xy)$neighs  }# find nn routine in splancs
    else {
     temp <- find.neighbor(xy,k=2);# find nn in spatialstats
     as.vector(temp[temp[,1] != temp[,2],2]);
# remove self neighbors, return id of neighbor
#  N.B. assumes find.neighbor returns list in order of points 1, 2,
     }
  }

