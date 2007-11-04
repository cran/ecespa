`ginv` <-
function(m) {

# compute a generalized inverse of matrix m

temp <- eigen(m, symmetric=T);# calculate eigenvalues and eigenvectors 

va <- temp$values;# extract the e-vals and e-vects
ve <- temp$vectors;

va <- ifelse((abs(va)<1e-9), 0, 1/va);
#  invert e-vals, but leave  0's alone

va2 <- 0*m;# matrix of zeros, same size as m
diag(va2) <- va# stuff into a square matrix

ve %*% va2 %*% t(ve);# and construct inverse 
}

