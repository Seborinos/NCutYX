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
# w.gaussian <- function(Z,Y,X,sigma2=1,w=Ps){
#   n <- dim(X)[1]
#   q <- dim(Z)[2]
#   p <- dim(Y)[2]
#   r <- dim(X)[2]
#   Dzy <- matrix(0,q,p)
#   Dyx <- matrix(0,p,r)
#   D <- matrix(0,q+p+r,q+p+r)
#   for (i in 1:q){
#     for (j in 1:p){
#       Dzy[i,j] <- exp((-1)*sum(Ps*(Z[,i]-Y[,j])^2))/sigma2)#This might be slow, do C++ instead
#       Dzy[j,i] <- Dzy[i,j]
#     }
#   }
#
#   for (i in 1:p){
#     for (j in 1:r){
#       Dyx[i,j] <- exp((-1)*sum(Ps*(Y[,i]-X[,j])^2))/sigma2)#This might be slow, do C++ instead
#       Dyx[j,i] <- Dyx[i,j]
#     }
#   }
#
#   D[1:q,(q+1):(q+p)] <- Dzy
#   D[(q+1):(q+p),1:q,] <- t(Dzy)
#   D[(q+1):(q+p),(q+p+1):(q+p+r)] <- Dyx
#   D[(q+p+1):(q+p+r),(q+1):(q+p)] <- t(Dyx)
#   return(D)
# }
