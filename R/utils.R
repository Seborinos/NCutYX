
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

l2n <- function(vec){
  return(sqrt(sum(vec^2)))
}

Prob <- function(argu1, argu2, L, b)
{
  T1 <- L*log(b+1)
  return(exp(-(argu2-argu1)/T1))
}

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

Kcs <- function(x){
  n <- length(x)
  K <- length(unique(x))
  Kcs <- matrix(0,n,K)
  for(i in 1:K){
    Kcs[which(x==i),i] <- 1
  }
  return(Kcs)
}

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
