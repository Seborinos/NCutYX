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

