`mNNinfo` <-
function(xy, label, nnid=NULL, splancs=TRUE) {

  #  data is xy coordinates in xy, labels in label

  # output is a list with ON = obs. contingency table, 
  #   EN = expected counts
  #   VarN = variance-covariance matrix of counts

  if (is.null(nnid)) {
    nnid <- NNid(xy,splancs);# get the nearest neighbor id
    }# if not already provided

  n <- table(label)
  N <- sum(n);

  l <- names(n);# names of the types of points
  k <- length(l);

  R <- sum((1:N) == nnid[nnid]);# number of refl NN
  Q1 <- table(c(1:6,table(nnid)))
  Q <- 2*Q1[2] + 6*Q1[3] + 12*Q1[4] + 20*Q1[5] + 30*Q1[6] - 70;
    # get # shared neighbors,  extra bit (1:6) is to ensure
     # that outer table() always generates a length 6 vector
    #  correct for the extra by subtracting (2+6+12+20+30)=70

  ON <- matrix(0,nrow=k, ncol=k); 

# loop over all pairs of types to get Obs count
  for (i in 1:k) { 
   for (j in 1:k) {
    ON[i,j] <- sum((label==l[i]) & (label[nnid] ==l[j]));
    }
   }

# use second function to calculate E N and Var N
  temp <- mNNinfo2(n, R, Q);
rownames(ON)<-colnames(ON)<- rownames(temp$EN)<- colnames(temp$EN)<-l


# then return all the bits
  list(ON=ON, EN=temp$EN, VarN=temp$VarN, R=R, Q=Q);
}

