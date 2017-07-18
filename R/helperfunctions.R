sample.cluster <- function(x){
  K <- length(x)
  J <- matrix(0,K,1)
  s1 <- sample(K,1,prob=x)
  J[s1,] <- 1
  return(J)
}