sample.cluster <- function(x){
  K <- length(x)
  J <- matrix(0,K,1)
  s1 <- sample(K,1,prob=x)
  J[s1,] <- 1
  return(J)
}

#Also maybe you need to scale the weights so that they equal 1
w.gauss <- function(Z,Y,X,sigma=1){
  n <- dim(X)[1]
  q <- dim(Z)[2]
  p <- dim(Y)[2]
  r <- dim(X)[2]
  m <- q + p + r
  D <- matrix(0,m,m)
  A <- exp((-sigma)*as.matrix(dist(t(cbind(Z,Y,X)),
                                   method='euclidean',
                                   diag=T,upper=T)))

  D[1:q,(q+1):(q+p)] <- A[1:q,(q+1):(q+p)]
  D[(q+1):(q+p),1:q] <- t(D[1:q,(q+1):(q+p)])
  D[(q+1):(q+p),(q+p+1):(q+p+r)] <- A[(q+1):(q+p),(q+p+1):(q+p+r)]
  D[(q+p+1):(q+p+r),(q+1):(q+p)] <- t(D[(q+1):(q+p),(q+p+1):(q+p+r)] )
  return(D)
}

#rewrite the function below in C++
#Also maybe you need to scale the weights so that they equal 1
w.cor <- function(Z,
                  Y,
                  X){
  n <- dim(X)[1]
  q <- dim(Z)[2]
  p <- dim(Y)[2]
  r <- dim(X)[2]
  m <- q + p + r
  D <- matrix(0,m,m)


  D[1:q,(q+1):(q+p)] <- abs(cor(Z,Y))
  D[(q+1):(q+p),1:q] <- t(D[1:q,(q+1):(q+p)])
  D[(q+1):(q+p),(q+p+1):(q+p+r)] <- abs(cor(Y,X))
  D[(q+p+1):(q+p+r),(q+1):(q+p)] <- t(D[(q+1):(q+p),(q+p+1):(q+p+r)] )
  return(D)
}

#' This function standardize a vector
#' @return The output is a standardized vector.
#' @param vec is a vector.
l2n <- function(vec){
  return(sqrt(sum(vec^2)))
}

#' This function calculates the probability to accept the updated result
#' when updated value of the objectove function is smaller than the old one.
#' Note that argu2 must be greater than argu1.
#' @return The output is the proobability to accept the updated result.
#' @param argu1 is the updates value of the objective function
#' @param argu2 is the old value of the objective function
#' @param b is the iteration times.
#' @param L is the temperature parameter.
Prob <- function(argu1, argu2, L, b)
{
  T1 <- L*log(b+1)
  return(exp(-(argu2-argu1)/T1))
}

#' This function calculates the weighed datasets for the AWNCut method.
#' @return A list with where WX is the weighed dataset X and WZ is the weighed dataset Z.
#' @param X is a n x p1 matrix
#' @param Z is a n x p2 matrix
#' @param ws is a vector of weights for both X and Z datasets.
AWNcut.W <- function(X, Z, ws){
  p1 <- ncol(X)
  p2 <- ncol(Z)
  #Calculate matrix W
  DistX <- as.matrix(dist(t(t(X)*ws[1:p1]), upper=T, diag=T))
  diag(DistX) <- 1
  WX <- DistX^(-1)
  DistZ <- as.matrix(dist(t(t(Z)*ws[(1+p1):(p1+p2)]), upper=T, diag=T))
  diag(DistZ) <- 1
  WZ <- DistZ^(-1)
  return(list(WX, WZ))
}

#' This function calculate the value of the objective function for the AWNCut method.
#' @return A list with where OP1 is the value of the NCut measure in X,
#' OP2 is the value of the NCut measure in Z,
#' TOP is the sum of the NCut measure in both X and Z,
#' Cor.X is a vector of the average correlation for X,
#' Cor.Z is a vector of the average correlation for Z and
#' Cor.perfeature is a combination of the average correlation for X and Z.
#' @param X is a n x p1 matrix.
#' @param Z is a n x p2 matrix.
#' @param WX is the weighed dataset X.
#' @param WZ is the weighed dataset Z.
#' @param Cs is clustering result.
#' @param tau is tuning parameter in the objective function.
AWNcut.OP <- function(X, Z, WX, WZ, Cs, tau)
{
  n <- nrow(X)
  p1 <- ncol(X)
  p2 <- ncol(Z)
  #Calculate cuts and cutvols for each dataset
  cutvolX <- NULL
  cutX <- NULL
  for(i in 1:ncol(Cs)){
    T1 <- WX[which(Cs[,i]==1),which(Cs[,i]==1)]
    T2 <- WX[which(Cs[,i]==1),which(Cs[,1]==0)]
    cutvolX <- c(cutvolX, sum(T1[upper.tri(T1)]))
    cutX <- c(cutX,sum(T2))
  }
  OP1 <- cutvolX/cutX

  cutvolZ <- NULL
  cutZ <- NULL
  for(i in 1:ncol(Cs)){
    T1 <- WZ[which(Cs[,i]==1),which(Cs[,i]==1)]
    T2 <- WZ[which(Cs[,i]==1),which(Cs[,1]==0)]
    cutvolZ <- c(cutvolZ, sum(T1[upper.tri(T1)]))
    cutZ <- c(cutZ,sum(T2))
  }
  OP2 <- cutvolZ/cutZ
  #Calculate average correlations per feature
  Cor.X <- 0
  Cor.Z <- 0
  for(i in 1:ncol(Cs)){
    T <- cor(X[which(Cs[,i]==1),],Z[which(Cs[,i]==1),])
    T[is.na(T)] <- 0
    Cor.X <- Cor.X+apply(abs(T),1,mean)
    Cor.Z <- Cor.Z+apply(abs(T),2,mean)
  }
  return(list(OP1            = sum(OP1),
              OP2            = sum(OP2),
              TOP            = sum(OP1) + tau*sum(OP2),
              Cor.X          = Cor.X,
              Cor.Z          = Cor.Z,
              Cor.perfeature = c(Cor.X,Cor.Z)))
}

#' This function updates the clustering result of the AWNCut method.
#' @return Cs is the updated clustering result
#' @param WX is the weighed dataset X.
#' @param WZ is the weighed dataset Z.
#' @param K is the number of clusters.
#' @param Cs is the clustering result from the las iteration.
AWNcut.UpdateCs <- function(WX, WZ, K, Cs)
{
  P <- NULL
  Cs.rec <- NULL
  for(i in 1:ncol(Cs)){
    # Calculate cutvol/cut for each dataset and record the pair
    # that has the largest Euclidean distance in each dataset.
    T1 <- WX[which(Cs[,i]==1),which(Cs[,i]==1)]
    T2 <- WX[which(Cs[,i]==1),which(Cs[,i]==0)]
    T3 <- WZ[which(Cs[,i]==1),which(Cs[,i]==1)]
    T4 <- WZ[which(Cs[,i]==1),which(Cs[,i]==0)]
    P <- c(P, sum(T1[upper.tri(T1)])/sum(T2)+sum(T3[upper.tri(T3)])/sum(T4))
    Cs.rec <- rbind(Cs.rec,cbind(which(T1==min(T1),arr.ind=T)[1,],which(T3==min(T3),arr.ind=T)[1,]))
  }
  P <- P/sum(P)
  #Select K(+) and K(-) according to cutvol/cut, K(+) is the first element in S and K(-) is the last element in S
  S <- sample(K,K,replace=F,prob=P)
  #Select the observation that needs to be rearranged and record its location in the clustering result of the last iteration
  Ob.rec <- sample(Cs.rec[S[1],],1)
  ob.rec <- which(Cs[,S[1]]==1)[Ob.rec]

  Cs.new <- Cs
  Cs.new[ob.rec,S[1]] <- 0
  Cs.new[ob.rec,S[K]] <- 1
  #To make sure each cluster has at least two observations
  t <- sum(apply(Cs.new,2,max))
  l <- min(apply(Cs.new,2,function(x){return(length(which(x==1)))}))
  if((t==K)&&(l>=2)){
    return(Cs.new)
  }else{
    return(Cs)
  }
}

#' This function updates the feature selection result of the AWNCut method.
#' @return The output is a vector of the standardized updated weights.
#' @param X is an n x p1 matrix.
#' @param Z is an n x p2 matrix.
#' @param K is the number of clusters.
#' @param WX is the weighed dataset X.
#' @param WZ is the weighed dataset Z.
#' @param b is the iteration times.
#' @param Cs is the old clustering result.
#' @param ws is the old feature selection result.
#' @param tau is a vector of tuning parameter tau in the objective function.
AWNcut.UpdateWs <- function(X, Z, K, WX, WZ, b, Cs, ws, tau)
{
  n <- nrow(X)
  p1 <- ncol(X)
  p2 <- ncol(Z)
  #Compute the pace of each iteration according to b
  pace.x <- 1/(10*ceiling(b/10)*sqrt(p1))
  pace.z <- 1/(10*ceiling(b/10)*sqrt(p2))
  pace1 <- rep(-1, p1)
  pace2 <- rep(-1, p2)
  #Calculate the average correlation of each feature
  temp <- AWNcut.OP(X, Z, WX, WZ, Cs, tau)$Cor.perfeature
  #Use Kmeans method to cluster the average weight into 2 clusters in each dataset
  temp1 <- kmeans(temp[1:p1],2)
  temp2 <- kmeans(temp[(1+p1):(p1+p2)],2)
  #For features in the cluster that has a higher level of average correlation, update their weight by adding the pace
  pace1[which(temp1$cluster==which.max(temp1$center))] <- 1
  pace2[which(temp2$cluster==which.max(temp2$center))] <- 1
  pace <- c(pace1,pace2)
  ws.new <- ws+c(pace*c(rep(pace.x,p1),rep(pace.z,p2)))
  #To make sure that each weight is nonngative
  for(i in 1:length(ws.new)){
    ws.new[i] <- max(ws.new[i],0)
  }
  return(c(ws.new[1:p1]/l2n(ws.new[1:p1]),ws.new[(p1+1):(p1+p2)]/l2n(ws.new[(p1+1):(p1+p2)])))
}

#' This function outputs the value of DBI for a clustering result.
#' @return The output is the value of DBI for a clustering result.
#' @param X is a n x p matrix with n observations and p variables.
#' @param K is the number of clusters.
#' @param Cs is a n x K matrix containing the clustering result of X. The entries
#'        in Cs is either 1 or 0. If the ith observation is in the kth cluster, then Cs(i,k)=1,
#'        otherwise, it equals to 0.
#' @param ws is a vector of the weights for the p variables. the length of ws equals p.
DBI <- function(X, K, Cs, ws){
  X <- scale(X)
  p1 <- ncol(X)
  w1 <- ws[1:p1]/l2n(ws[1:p1])
  WX <- as.matrix(dist(t(t(X)*w1), upper=T, diag=T))

  #Calculate both cut and cuvol, cutvol is the diagonal of the matrix Cut, cut is the reat of the entries
  Cut <- matrix(0,K,K)
  for(i in 1:K){
    for(j in 1:K){
      T <- WX[which(Cs[,i]==1),which(Cs[,j]==1)]
      Cut[i,j] <- sum(T)
    }
  }
  cutvol <- diag(Cut)/2
  DBI <- NULL
  for(i in 1:K){
    t <- NULL
    for(j in (1:K)[-i]){
      t <- c(t,Cut[i,j]/(cutvol[i]+cutvol[j]))
    }
    DBI <- c(DBI, max(t))
  }
  return(mean(DBI[DBI<Inf])) #When a clustering result contains two clusters that has only one observation, the value of DBI will be Inf.
}

#' This function transfers a clustering result into a matrix format from a vector.
#' @return A clustering result in the matrix format.
#' @param x is a clustering result in the vector format.
Kcs <- function(x){
  n <- length(x)
  K <- length(unique(x))
  Kcs <- matrix(0,n,K)
  for(i in 1:K){
    Kcs[which(x==i),i] <- 1
  }
  return(Kcs)
}

#' This function calculates the stablility of the simulation result
#' @return The output if the stability of the simulation result
#' @param x is a 3 dimensional array.
#'          The first dimension equals to sample size.
#'          The second dimension equals to number og clusters.
#'          The third dimension equals to the replication times.
Stability <- function(x){
  N <- dim(x)[3]
  n <- dim(x)[1]
  sta <- 0
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      sta <- sta+sum(abs(x[,,i]-x[,,j]))/(n^2)
    }
  }
  return(sta/choose(N,2))
}
