sample.cluster <- function(x){
  K <- length(x)
  J <- matrix(0,K,1)
  s1 <- sample(K,1,prob=x)
  J[s1,] <- 1
  return(J)
}

neovector <- function(Probs){
  K <- length(Probs)
  return(rbinom(K,1,prob=Probs))
}


#rewrite the function below in C++
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


OP <- function(X,
               Z,
               WX,
               WZ,
               Cs){

  n  <- nrow(X)
  p1 <- ncol(X)
  p2 <- ncol(Z)
  # Cut and cutvol
  cutvolX <- NULL
  cutX    <- NULL
  for(i in 1:ncol(Cs)){
    T1      <- WX[which(Cs[,i]==1),which(Cs[,i]==1)]
    T2      <- WX[which(Cs[,i]==1),which(Cs[,1]==0)]
    cutvolX <- c(cutvolX, sum(T1[upper.tri(T1)]))
    cutX    <- c(cutX,sum(T2))
  }
  OP1 <- cutvolX/cutX

  cutvolZ <- NULL
  cutZ    <- NULL
  for(i in 1:ncol(Cs)){
    T1      <- WZ[which(Cs[,i]==1),which(Cs[,i]==1)]
    T2      <- WZ[which(Cs[,i]==1),which(Cs[,1]==0)]
    cutvolZ <- c(cutvolZ, sum(T1[upper.tri(T1)]))
    cutZ    <- c(cutZ,sum(T2))
  }
  OP2 <- cutvolZ/cutZ
  # Correlations perfeature
  Cor.X <- 0
  Cor.Z <- 0
  for(i in 1:ncol(Cs)){
    T           <- cor(X[which(Cs[,i]==1),],Z[which(Cs[,i]==1),])
    T[is.na(T)] <- 0
    Cor.X       <- Cor.X+apply(abs(T),1,mean)
    Cor.Z       <- Cor.Z+apply(abs(T),2,mean)
  }
  return(list(OP1=sum(OP1),
              OP2=sum(OP2),
              TOP=sum(OP1)+sum(OP2),
              Cor.X=Cor.X, Cor.Z=Cor.Z,
              Cor.perfeature=c(Cor.X,Cor.Z)))
}


