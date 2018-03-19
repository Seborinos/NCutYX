#' Cluster the Columns of Y into K Groups Using the NCut Graph Measure.
#'
#' Builds a similarity matrix for the columns of Y and clusters them into
#' K groups based on the NCut graph measure. Correlation, Euclidean and Gaussian distances can be used
#' to construct the similarity matrix.
#' @param Y is a n x p matrix of p variables and n observations. The p columns of
#' Y will be clustered into K groups using NCut.
#' @param K is the number of clusters.
#' @param B is the number of iterations.
#' @param N is the number of samples per iterations.
#' @param dist is the type of distance metric for the construction of the similarity matrix.
#' Options are 'gaussian', 'euclidean' and 'correlation', the latter being the default.
#' @param scale equals TRUE if data Y is to be scaled with mean 0 and variance 1.
#' @param q is the quantile used for the top results at each iterations.
#' @param sigma is the bandwidth parameter when the dist metric chosen is 'gaussian' (default=0.1).
#' @return A list with the following components:
#' \describe{
#' \item{quantile}{a vector of length \code{N} which contains the quantiles
#' \code{q} at each iteration of the optimization algorithm.}
#' \item{cluster}{a matrix representing the clustering result of dimension \code{p} times
#' \code{K}, where \code{p} is the number of columns of \code{Y}.}
#' \item{ncut}{the NCut measure for the cluster result.}
#' }
#' @details
#' The algorithm minimizes the NCut through the cross entropy method.
#' The edges of the graph correspond to the entries of a similarity matrix constructed based on a
#' correlation, euclidean or gaussian distance metric.
#' The clusters correspond to partitions that minimize this NCut objective function.
#' @author Sebastian Jose Teran Hidalgo. Maintainer: Sebastian Jose Teran Hidalgo
#' \url{sebastianteranhidalgo@gmail.com}.
#' @references Von Luxburg, Ulrike. "A tutorial on spectral clustering."
#' Statistics and computing 17.4 (2007): 395-416.
#'
#' Kroese, D. P., Rubinstein, R. Y., Cohen, I., Porotsky, S., & Taimre, T. (2013).
#' "Cross-entropy method."
#' In Encyclopedia of Operations Research and Management Science (pp. 326-333). Springer US.
#' @examples
#' # This sets up the initial parameters for the simulation.
#' library(MASS)
#' n=100 # Sample size
#' B=30 # Number of iterations in the simulated annealing algorithm.
#' p=50 # Number of columns of Y.
#'
#' S=matrix(0.2,p,p)
#' S[1:(p/2),(p/2+1):p]=0
#' S[(p/2+1):p,1:(p/2)]=0
#' S=S-diag(diag(S))+diag(p)
#' mu=rep(0,p)
#'
#' W0=matrix(1,p,p)
#' W0[1:(p/2),1:(p/2)]=0
#' W0[(p/2+1):p,(p/2+1):p]=0
#' Denum=sum(W0)
#'
#' Y=mvrnorm(n, mu, S)
#' # NCut
#' Res=ncut(Y,
#' K=2,
#' B=30,
#' N=1000,
#' dist='correlation',
#' scale=TRUE,
#' q=0.2,
#' sigma=0.1)
#' Cx=Res[[2]]
#' f11=matrix(Cx[,1],p,1)
#' f12=matrix(Cx[,2],p,1)
#'
#' errorL=sum((f11%*%t(f11))*W0)/Denum+sum((f12%*%t(f12))*W0)/Denum
#' # This is the true error of the clustering solution.
#' errorL
#' @export
ncut <- function(Y,
                 K     = 2,
                 B     = 30,
                 N     = 500,
                 dist  = 'correlation',
                 scale = TRUE,
                 q     = 0.1,
                 sigma = 1){
  # This creates the weight matrix W
  Res <- list()
  quantiles <- vector(mode="numeric", length=B)
  if(scale==T){
    Y <- scale(Y)
  }
  p <- dim(Y)[2]
  if(dist=='euclidean'){
    Wyy <- as.matrix(stats::dist(t(Y),diag=T,upper=T)) + diag(p)
    Wyy <- Wyy^(-1)
  }else if(dist=='gaussian'){
    Wyy <- exp((-1)*as.matrix(stats::dist(t(Y),diag=T,upper=T))/sigma)
  }else if(dist=='correlation'){
    Wyy <- abs(stats::cor(Y))
  }else{
    print('Distance Error')
  }
  # vector with probabilities of mus being 0 or not
  Ps <- matrix(1/K,p,K)
  q0 <- as.integer(q*N)
  p0 <- 1/K
  return(ncutcem(Wyy,p,K,N,B,q0,p0))
}

#' Cluster the Columns of Y into K Groups with the Help of External Features X.
#'
#' This function will output K clusters of the columns of Y using the help of
#' X.
#' @param Y is a n x p matrix of p variables and n observations. The columns of
#' Y will be clustered into K groups.
#' @param X is a n x q matrix of q variables and n observations.
#' @param K is the number of clusters.
#' @param sampling if 'equal' then the sampling probabilities is the same during
#' the simulated annealing algorithm, if 'size' the probabilites are proportional
#' the the sizes of the clusters in the current iterations.
#' @param B is the number of iterations in the simulated annealing algorithm.
#' @param L is the temperature coefficient in the simulated annealing algorithm.
#' @param alpha is the coefficient of the elastic net penalty.
#' @param nlambdas is the number of tuning parameters in the elastic net.
#' @param ncv is the number of cross-validations in the elastic net.
#' @param dist is the type of distance metric for the construction of the similarity matrix.
#' Options are 'gaussian', 'euclidean' and 'correlation', the latter being the default.
#' @param sigma is the parameter for the gaussian kernel distance which is ignored if 'gaussian' is
#' not chosen as distance measure.
#' @return A list with the final value of the objective function,
#' the clusters and the lambda penalty chosen through cross-validation.
#' @details
#' The algorithm minimizes a modified version of NCut through simulated annealing.
#' The modified NCut uses in the numerator the similarity matrix of the original data \code{Y}
#' and the denominator uses the similarity matrix of the prediction of \code{Y} using \code{X}.
#' The clusters correspond to partitions that minimize this objective function.
#' The external information of \code{X} is incorporated by using elastic net to predict \code{Y}.
#' @author Sebastian Jose Teran Hidalgo and Shuangge Ma. Maintainer: Sebastian Jose Teran Hidalgo.
#' \url{sebastianteranhidalgo@gmail.com}.
#' @references Hidalgo, Sebastian J. Teran, Mengyun Wu, and Shuangge Ma.
#' Assisted clustering of gene expression data using ANCut. BMC genomics 18.1 (2017): 623.
#' @return A list with the following components:
#' \describe{
#' \item{loss}{a vector of length \code{N} which contains the loss
#' at each iteration of the simulated annealing algorithm.}
#' \item{cluster}{a matrix representing the clustering result of dimension \code{p} times
#' \code{K}, where \code{p} is the number of columns of \code{Y}.}
#' \item{lambda.min}{is the optimal lambda chosen through cross-validation for the elastic net for
#' predicting \code{Y} with \code{Y}.}
#' }
#' @examples
#' #This sets up the initial parameters for the simulation.
#' library(MASS)#for mvrnorm
#' library(fields)
#' n=30 #Sample size
#' B=50 #Number of iterations in the simulated annealing algorithm.
#' L=10000 #Temperature coefficient.
#' p=50 #Number of columns of Y.
#' q=p #Number of columns of X.
#' h1=0.15
#' h2=0.25
#'
#' S=matrix(0.2,q,q)
#' S[1:(q/2),(q/2+1):q]=0
#' S[(q/2+1):q,1:(q/2)]=0
#' S=S-diag(diag(S))+diag(q)
#'
#' mu=rep(0,q)
#'
#' W0=matrix(1,p,p)
#' W0[1:(p/2),1:(p/2)]=0
#' W0[(p/2+1):p,(p/2+1):p]=0
#' Denum=sum(W0)
#'
#' B2=matrix(0,q,p)
#' for (i in 1:(p/2)){
#'    B2[1:(q/2),i]=runif(q/2,h1,h2)
#'    in1=sample.int(q/2,6)
#'    B2[-in1,i]=0
#' }
#'
#' for (i in (p/2+1):p){
#'    B2[(q/2+1):q,i]=runif(q/2,h1,h2)
#'    in2=sample(seq(q/2+1,q),6)
#'    B2[-in2,i]=0
#' }
#'
#' X=mvrnorm(n, mu, S)
#' Z=X%*%B2
#' Y=Z+matrix(rnorm(n*p,0,1),n,p)
#' #Our method
#' Res=ancut(Y=Y,X=X,B=B,L=L,alpha=0,ncv=3)
#' Cx=Res[[2]]
#' f11=matrix(Cx[,1],p,1)
#' f12=matrix(Cx[,2],p,1)
#'
#' errorL=sum((f11%*%t(f11))*W0)/Denum+sum((f12%*%t(f12))*W0)/Denum
#' #This is the true error of the clustering solution.
#' errorL
#'
#' par(mfrow=c(1,2))
#' #Below is a plot of the simulated annealing path.
#' plot(Res[[1]],type='l',ylab='')
#' #Cluster found by ANCut
#' image.plot(Cx)
#' @export
ancut <- function(Y,
                  X,
                  K        = 2,
                  B        = 3000,
                  L        = 1000,
                  alpha    = 0.5,
                  nlambdas = 100,
                  sampling = 'equal',
                  ncv      = 5,
                  dist     = 'correlation',
                  sigma    = 0.1){
  # This creates the weight matrix W
  X <- scale(X)
  Y <- scale(Y)
  p <- dim(Y)[2]

  p <- dim(Y)[2]
  if(dist=='euclidean'){
    Wyy <- as.matrix(stats::dist(t(Y),diag=T,upper=T)) + diag(p)
    Wyy <- Wyy^(-1)
  }else if(dist=='gaussian'){
    Wyy <- exp((-1)*as.matrix(stats::dist(t(Y),diag=T,upper=T))/sigma)
  }else if(dist=='correlation'){
    Wyy <- abs(stats::cor(Y))
  }else{
    print('Distance Error')
  }

  #modelling the relationship between Y and X
  cv.m1 <- glmnet::cv.glmnet(X,
                             Y,
                             family    = c("mgaussian"),
                             alpha     = alpha,
                             nfolds    = ncv,
                             nlambda   = nlambdas,
                             intercept = FALSE)

  m1 <- glmnet::glmnet(X,
                       Y,
                       family    = c("mgaussian"),
                       alpha     = alpha,
                       lambda    = cv.m1$lambda.min,
                       intercept = FALSE)

  Y2 <- predict(m1,newx=X)
  Y2 <- scale(Y2[ , ,1])
  if(dist=='euclidean'){
    Wxx <- as.matrix(stats::dist(t(Y2),diag=T,upper=T)) + diag(p)
    Wxx <- Wxx^(-1)
  }else if(dist=='gaussian'){
    Wxx <- exp((-1)*as.matrix(stats::dist(t(Y2),diag=T,upper=T))/sigma)
  }else if(dist=='correlation'){
    Wxx <- abs(stats::cor(Y2))
  }else{
    print('Distance Error')
  }

  #This creates a random starting point in the split in the algorithm for K clusters
  Cx <- matrix(0,p,K)

  for (i in 1:p){
    Cx[i,sample(K,1)] <- 1
  }

  #Now, calculate the number of indices in each group.
  Nx <- apply(Cx[,1:K],2,sum)

  #These matrices will keep track of the elements of the clusters while
  #doing simulated annealing.
  C2x  <- matrix(0,p,K+1)
  C2x  <- Cx
  J    <- NCutY3V1(Cx[ ,1:K],matrix(1,p,K)-Cx[ ,1:K],Wyy,Wxx)
  Test <- vector(mode="numeric", length=B)

  for (k in 1:B){
    ###Draw k(-) and k(+)with unequal probabilites.
    #This section needs to change dramatically for
    #the general case
    N <- sum(Nx)
    P <- Nx/N

    if(sampling=='equal'){
      s <- sample.int(K,K,replace=FALSE)
    }else if(sampling=='size'){
      s <- sample.int(K,K,replace=FALSE,prob=P)
    }

    ###Select a vertex from cluster s[1] with unequal probability
    #Calculating Unequal probabilites
    #Draw a coin to see whether we choose X or Y
    ax           <- which(Cx[ ,s[1]]==1)#which Xs belong to the cluster
    sx           <- sample(ax,1)
    C2x[sx,s[1]] <- 0
    C2x[sx,s[K]] <- 1

    #Now Step 3 in the algorithm
    J2 <- NCutY3V1(C2x[ ,1:K],matrix(1,p,K)-C2x[ ,1:K],Wyy,Wxx)

    if (J2>J){
      des <- rbinom(1,1,exp(-L*log(k+1)*(J2-J)))
      if (des==1){
        Cx <- C2x#Set-up the new clusters
        J  <- J2
        Nx <- apply(Cx[ ,1:K],2,sum)
      }else{
        C2x <- Cx
      }
    } else{
      Cx <- C2x
      J  <- J2
      Nx <- apply(Cx[ ,1:K],2,sum)
    }
    Test[k] <- J

  }
  Res            <- list()
  Res$loss       <- Test
  Res$cluster    <- Cx
  Res$lambda.min <- cv.m1$lambda.min
  return(Res)
}

#' MuNCut Clusters the Columns of Data from 3 Different Sources.
#'
#' It clusters the columns of Z,Y and X into K clusters by representing each data type as one network layer.
#' It represents the Z layer depending on Y, and the Y layer depending on X. Elastic net can be used before the clustering
#' procedure by using the predictions of Z and Y instead of the actual values to improve the cluster results.
#' This function will output K clusters of columns of Z, Y and X.
#' @param Z is a n x q matrix of q variables and n observations.
#' @param Y is a n x p matrix of p variables and n observations.
#' @param X is a n x r matrix of r variables and n observations.
#' @param K is the number of column clusters.
#' @param B is the number of iterations in the simulated annealing algorithm.
#' @param L is the temperature coefficient in the simulated annealing algorithm.
#' @param alpha is the tuning parameter in the elastic net penalty, only used when model=T.
#' @param ncv is the number of cross-validations used to choose the tuning parameter lambda in the
#' elastic net penalty, only used when model=T.
#' @param nlambdas number of tuning parameters lambda used during cross-validation, only when model=T.
#' @param scale when TRUE the Z, Y and X are scaled with mean 0 and standard deviation equal 1.
#' @param model when TRUE the the relationship between Z and Y, and between Y and X are modeled
#' with the elastic net. The predictions of Z, and Y from the models are used in the clustering algorithm.
#' @param gamma is the tuning parameter of the clustering penalty. Larger values give more importance to
#' within layer effects and less to across layer effects.
#' @param sampling if 'equal' then the sampling distribution is discrete uniform over the
#' number of clusters, if 'size' the probabilities are inversely proportional to the size
#' of each cluster.
#' @param dist is the type of distance measure use in the similarity matrix.
#' Options are 'gaussian' and 'correlation', with 'gaussian' being the default.
#' @param sigma is the bandwidth parameter when the dist metric chosen is gaussian.
#' @details
#' The algorithm minimizes a modified version of NCut through simulated annealing.
#' The clusters correspond to partitions that minimize this objective function.
#' The external information of X is incorporated by using ridge regression to predict Y.
#' @references  Sebastian J. Teran Hidalgo and Shuangge Ma.
#' Clustering Multilayer Omics Data using MuNCut. (Revise and resubmit.)
#' @examples
#' library(NCutYX)
#' library(MASS)
#' library(fields) #for image.plot
#'
#' #parameters#
#' set.seed(777)
#' n=50
#' p=50
#' h=0.5
#' rho=0.5
#'
#' W0=matrix(1,p,p)
#' W0[1:(p/5),1:(p/5)]=0
#' W0[(p/5+1):(3*p/5),(p/5+1):(3*p/5)]=0
#' W0[(3*p/5+1):(4*p/5),(3*p/5+1):(4*p/5)]=0
#' W0[(4*p/5+1):p,(4*p/5+1):p]=0
#' W0=cbind(W0,W0,W0)
#' W0=rbind(W0,W0,W0)
#'
#' Y=matrix(0,n,p)
#' Z=matrix(0,n,p)
#' Sigma=matrix(rho,p,p)
#' Sigma[1:(p/5),1:(p/5)]=2*rho
#' Sigma[(p/5+1):(3*p/5),(p/5+1):(3*p/5)]=2*rho
#' Sigma[(3*p/5+1):(4*p/5),(3*p/5+1):(4*p/5)]=2*rho
#' Sigma=Sigma-diag(diag(Sigma))
#' Sigma=Sigma+diag(p)
#'
#' X=mvrnorm(n,rep(0,p),Sigma)
#' B1=matrix(0,p,p)
#' B2=matrix(0,p,p)
#'
#' B1[1:(p/5),1:(p/5)]=runif((p/5)^2,h/2,h)*rbinom((p/5)^2,1,0.2)
#' B1[(p/5+1):(3*p/5),(p/5+1):(3*p/5)]=runif((2*p/5)^2,h/2,h)*rbinom((2*p/5)^2,1,0.2)
#' B1[(3*p/5+1):(4*p/5),(3*p/5+1):(4*p/5)]=runif((p/5)^2,h/2,h)*rbinom((p/5)^2,1,0.2)
#'
#' B2[1:(p/5),1:(p/5)]=runif((p/5)^2,h/2,h)*rbinom((p/5)^2,1,0.2)
#' B2[(p/5+1):(3*p/5),(p/5+1):(3*p/5)]=runif((2*p/5)^2,h/2,h)*rbinom((2*p/5)^2,1,0.2)
#' B2[(3*p/5+1):(4*p/5),(3*p/5+1):(4*p/5)]=runif((p/5)^2,h/2,h)*rbinom((p/5)^2,1,0.2)
#'
#' Y=X%*%B1+matrix(rnorm(n*p,0,0.5),n,p)
#' Y2=X%*%B1
#'
#' Z=Y%*%B2+matrix(rnorm(n*p,0,0.5),n,p)
#' Z2=Y%*%B2
#'
#' #Computing our method
#' clust <- muncut(Z,
#'                 Y,
#'                 X,
#'                 K        = 4,
#'                 B        = 10000,
#'                 L        = 500,
#'                 sampling = 'size',
#'                 alpha    = 0.5,
#'                 ncv      = 3,
#'                 nlambdas = 20,
#'                 sigma    = 10,
#'                 scale    = TRUE,
#'                 model    = FALSE,
#'                 gamma    = 0.1)
#'
#' A <- clust[[2]][,1]%*%t(clust[[2]][,1])+
#'      clust[[2]][,2]%*%t(clust[[2]][,2])+
#'      clust[[2]][,3]%*%t(clust[[2]][,3])+
#'      clust[[2]][,4]%*%t(clust[[2]][,4])
#'
#' errorK=sum(A*W0)/(3*p)^2
#' errorK
#' plot(clust[[1]],type='l')
#' image.plot(A)
#' @export
muncut <- function(Z,
                   Y,
                   X,
                   K        = 2,
                   B        = 3000,
                   L        = 1000,
                   alpha    = 0.5,
                   ncv      = 3,
                   nlambdas = 100,
                   scale    = FALSE,
                   model    = FALSE,
                   gamma    = 0.5,
                   sampling = 'equal',
                   dist     = 'gaussian',
                   sigma    = 0.1){
  # Beginning of the function
  if (model==T){
    if (scale==T){
      Z <- scale(Z)
      Y <- scale(Y)
      X <- scale(X)
      q <- dim(Z)[2]
      p <- dim(Y)[2]
      r <- dim(X)[2]
      m <- q + p + r
      #Elastic net to predict Y with X
      cv.m1 <- glmnet::cv.glmnet(X,
                                 Y,
                                 family=c("mgaussian"),
                                 alpha=alpha,
                                 nfolds=ncv,
                                 nlambda=nlambdas,
                                 intercept=FALSE)

      m1    <- glmnet::glmnet(X,
                              Y,
                              family=c("mgaussian"),
                              alpha=alpha,
                              lambda=cv.m1$lambda.min,
                              intercept=FALSE)

      Y2 <- predict(m1,newx=X)
      Y2 <- scale(Y2[ , ,1])
      #Elastic net to predict Z with Y
      cv.m2 <- glmnet::cv.glmnet(Y,
                                 Z,
                                 family    = c("mgaussian"),
                                 alpha     = alpha,
                                 nfolds    = ncv,
                                 nlambda   = nlambdas,
                                 intercept = FALSE)

      m2    <- glmnet::glmnet(Y,
                              Z,
                              family    = c("mgaussian"),
                              alpha     = alpha,
                              lambda    = cv.m1$lambda.min,
                              intercept = FALSE)

      Z2 <- predict(m2,newx=Y)
      Z2 <- scale(Z2[ , ,1])
      if (dist=='euclidean'){
        #Distance matrix for the predicted variables
        ZYX2  <- cbind(Z2,Y2,X)
        Wzyx2 <- as.matrix(stats::dist(t(ZYX2),diag=T,upper=T)) + diag(m)
        Wzyx2 <- Wzyx2^(-1)
        #Z's distance matrix
        Wz    <- as.matrix(stats::dist(t(Z),diag=T,upper=T)) + diag(q)
        Wz    <- Wz^(-1)
        #Y's distance matrix
        Wy    <- as.matrix(stats::dist(t(Y),diag=T,upper=T)) + diag(p)
        Wy    <- Wy^(-1)
        #X's distance matrix
        Wx    <- as.matrix(stats::dist(t(X),diag=T,upper=T)) + diag(r)
        Wx    <- Wx^(-1)
        #Matrix without diagonal entries
        Izyx2                      <- Wzyx2
        Izyx2[1:q,1:q]             <- 0
        Izyx2[(1:p+q),(1:p+q)]     <- 0
        Izyx2[(1:r+p+q),(1:r+p+q)] <- 0
      }else if(dist=='gaussian'){
        #Distance matrix for the predicted variables
        ZYX2  <- cbind(Z2,Y2,X)
        Wzyx2 <- exp((-1)*as.matrix(stats::dist(t(ZYX2),diag=T,upper=T))/sigma)
        #Z's distance matrix
        Wz    <- exp((-1)*as.matrix(stats::dist(t(Z),diag=T,upper=T))/sigma)
        #Y's distance matrix
        Wy    <- exp((-1)*as.matrix(stats::dist(t(Y),diag=T,upper=T))/sigma)
        #X's distance matrix
        Wx    <- exp((-1)*as.matrix(stats::dist(t(X),diag=T,upper=T))/sigma)
        #Matrix without diagonal entries
        Izyx2                      <- Wzyx2
        Izyx2[1:q,1:q]             <- 0
        Izyx2[(1:p+q),(1:p+q)]     <- 0
        Izyx2[(1:r+p+q),(1:r+p+q)] <- 0
      } else{
        print('Distance Error')
      }
    }else{
      q <- dim(Z)[2]
      p <- dim(Y)[2]
      r <- dim(X)[2]
      m <- q + p + r
      #Elastic net to predict Y with X
      cv.m1 <- glmnet::cv.glmnet(X,
                                 Y,
                                 family    = c("mgaussian"),
                                 alpha     = alpha,
                                 nfolds    = ncv,
                                 nlambda   = nlambdas,
                                 intercept = FALSE)
      m1    <- glmnet::glmnet(X,
                              Y,
                              family    = c("mgaussian"),
                              alpha     = alpha,
                              lambda    = cv.m1$lambda.min,
                              intercept = FALSE)
      Y2 <- predict(m1,newx=X)
      Y2 <- Y2[ , ,1]
      #Elastic net to predict Z with Y
      cv.m2 <- glmnet::cv.glmnet(Y,
                                 Z,
                                 family    = c("mgaussian"),
                                 alpha     = alpha,
                                 nfolds    = ncv,
                                 nlambda   = nlambdas,
                                 intercept = FALSE)
      m2    <- glmnet::glmnet(Y,
                              Z,
                              family    = c("mgaussian"),
                              alpha     = alpha,
                              lambda    = cv.m2$lambda.min,
                              intercept = FALSE)
      Z2 <- predict(m2,newx=Y)
      Z2 <- Z2[ , ,1]
      if (dist=='euclidean'){
        #Distance matrix for the predicted variables
        ZYX2  <- cbind(Z2,Y2,X)
        Wzyx2 <- as.matrix(stats::dist(t(ZYX2),diag=T,upper=T)) + diag(q+p+r)
        Wzyx2 <- Wzyx2^(-1)
        #Z's distance matrix
        Wz <- as.matrix(stats::dist(t(Z),diag=T,upper=T)) + diag(q)
        Wz <- Wz^(-1)
        #Y's distance matrix
        Wy <- as.matrix(stats::dist(t(Y),diag=T,upper=T)) + diag(p)
        Wy <- Wy^(-1)
        #X's distance matrix
        Wx <- as.matrix(stats::dist(t(X),diag=T,upper=T)) + diag(r)
        Wx <- Wx^(-1)
        #Matrix without diagonal entries
        Izyx2                      <- Wzyx2
        Izyx2[1:q,1:q]             <- 0
        Izyx2[(1:p+q),(1:p+q)]     <- 0
        Izyx2[(1:r+p+q),(1:r+p+q)] <- 0
      }else if(dist=='gaussian'){
        #Distance matrix for the predicted variables
        ZYX2  <- cbind(Z2,Y2,X)
        Wzyx2 <- exp((-1)*as.matrix(stats::dist(t(ZYX2),diag=T,upper=T))/sigma)
        #Z's distance matrix
        Wz <- exp((-1)*as.matrix(stats::dist(t(Z),diag=T,upper=T))/sigma)
        #Y's distance matrix
        Wy <- exp((-1)*as.matrix(stats::dist(t(Y),diag=T,upper=T))/sigma)
        #X's distance matrix
        Wx <- exp((-1)*as.matrix(stats::dist(t(X),diag=T,upper=T))/sigma)
        #Matrix without diagonal entries
        Izyx2                      <- Wzyx2
        Izyx2[1:q,1:q]             <- 0
        Izyx2[(1:p+q),(1:p+q)]     <- 0
        Izyx2[(1:r+p+q),(1:r+p+q)] <- 0
      } else{
        print('Distance Error')
      }
    }
  }else{
    if (scale==T){
      Z <- scale(Z)
      Y <- scale(Y)
      X <- scale(X)
      q <- dim(Z)[2]
      p <- dim(Y)[2]
      r <- dim(X)[2]
      if (dist=='euclidean'){
        #Distance matrix for the predicted variables
        ZYX2  <- cbind(Z,Y,X)
        Wzyx2 <- as.matrix(stats::dist(t(ZYX2),diag=T,upper=T)) + diag(q+p+r)
        Wzyx2 <- Wzyx2^(-1)
        #Z's distance matrix
        Wz <- as.matrix(stats::dist(t(Z),diag=T,upper=T)) + diag(q)
        Wz <- Wz^(-1)
        #Y's distance matrix
        Wy <- as.matrix(stats::dist(t(Y),diag=T,upper=T)) + diag(p)
        Wy <- Wy^(-1)
        #X's distance matrix
        Wx <- as.matrix(stats::dist(t(X),diag=T,upper=T)) + diag(r)
        Wx <- Wx^(-1)
        #Matrix without diagonal entries
        Izyx2                      <- Wzyx2
        Izyx2[1:q,1:q]             <- 0
        Izyx2[(1:p+q),(1:p+q)]     <- 0
        Izyx2[(1:r+p+q),(1:r+p+q)] <- 0
      }else if(dist=='gaussian'){
        #Distance matrix for the predicted variables
        ZYX2  <- cbind(Z,Y,X)
        Wzyx2 <- exp((-1)*as.matrix(dist(t(ZYX2),diag=T,upper=T))/sigma)
        #Z's distance matrix
        Wz <- exp((-1)*as.matrix(dist(t(Z),diag=T,upper=T))/sigma)
        #Y's distance matrix
        Wy <- exp((-1)*as.matrix(dist(t(Y),diag=T,upper=T))/sigma)
        #X's distance matrix
        Wx <- exp((-1)*as.matrix(dist(t(X),diag=T,upper=T))/sigma)
        #Matrix without diagonal entries
        Izyx2                      <- Wzyx2
        Izyx2[1:q,1:q]             <- 0
        Izyx2[(1:p+q),(1:p+q)]     <- 0
        Izyx2[(1:r+p+q),(1:r+p+q)] <- 0
      } else{
        print('Distance Error')
      }
    }else{
      q <- dim(Z)[2]
      p <- dim(Y)[2]
      r <- dim(X)[2]
      if (dist=='euclidean'){
        #Distance matrix for the predicted variables
        ZYX2  <- cbind(Z,Y,X)
        Wzyx2 <- as.matrix(dist(t(ZYX2),diag=T,upper=T)) + diag(q+p+r)
        Wzyx2 <- Wzyx2^(-1)
        #Z's distance matrix
        Wz <- as.matrix(dist(t(Z),diag=T,upper=T)) + diag(q)
        Wz <- Wz^(-1)
        #Y's distance matrix
        Wy <- as.matrix(dist(t(Y),diag=T,upper=T)) + diag(p)
        Wy <- Wy^(-1)
        #X's distance matrix
        Wx <- as.matrix(dist(t(X),diag=T,upper=T)) + diag(r)
        Wx <- Wx^(-1)
        #Matrix without diagonal entries
        Izyx2                      <- Wzyx2
        Izyx2[1:q,1:q]             <- 0
        Izyx2[(1:p+q),(1:p+q)]     <- 0
        Izyx2[(1:r+p+q),(1:r+p+q)] <- 0
      }else if(dist=='gaussian'){
        #Distance matrix for the predicted variables
        ZYX2  <- cbind(Z,Y,X)
        Wzyx2 <- exp((-1)*as.matrix(dist(t(ZYX2),diag=T,upper=T))/sigma)
        #Z's distance matrix
        Wz <- exp((-1)*as.matrix(dist(t(Z),diag=T,upper=T))/sigma)
        #Y's distance matrix
        Wy <- exp((-1)*as.matrix(dist(t(Y),diag=T,upper=T))/sigma)
        #X's distance matrix
        Wx <- exp((-1)*as.matrix(dist(t(X),diag=T,upper=T))/sigma)
        #Matrix without diagonal entries
        Izyx2                      <- Wzyx2
        Izyx2[1:q,1:q]             <- 0
        Izyx2[(1:p+q),(1:p+q)]     <- 0
        Izyx2[(1:r+p+q),(1:r+p+q)] <- 0
      } else{
        print('Distance Error')
      }
    }
  }

  #This creates a random starting point in the split in the algorithm for K clusters
  Cx <- matrix(0,q+p+r,K)
  #This is a matrix of only ones
  M1 <- matrix(1,q+p+r,K)
  #Below: force the result to have one member per type of data and cluster.
  #The code below makes sure each cluster gets at least Min elements
  #per data type
  Min   <- 2
  Check <- matrix(0,3,K)
  while(sum((Check<=Min))>0){
    for (i in 1:(q+p+r)){
      Cx[i,sample(K,1)] <- 1
    }
    Check[1, ] <- apply(Cx[1:q,1:K], 2, sum)
    Check[2, ] <- apply(Cx[(q+1):(q+p),1:K], 2, sum)
    Check[3, ] <- apply(Cx[(p+q+1):(q+p+r),1:K], 2, sum)
  }

  #Now, calculate the number of indices in each group.
  Nx <- apply(Cx[ ,1:K], 2, sum)

  #These matrices will keep track of the elements of the clusters while
  #doing simulated annealing.
  C2x <- matrix(0, p+q+r, K)
  C2x <- Cx
  J   <- NCut(Cx[,1:K], Izyx2) + gamma*(NCut(Cx[1:q,1:K], Wz) +
         NCut(Cx[(q+1:p),1:K], Wy) + NCut(Cx[(q+p+1:r),1:K], Wx))

  # J=NCutY3V1(Cx[,1:K],M1-Cx[,1:K],Izyx2,Izyx2)+
  #   gamma*(NCutY3V1(Cx[1:q,1:K],M1[1:q]-Cx[1:q,1:K],Wz,Wz)+
  #   NCutY3V1(Cx[(q+1:p),1:K],M1[1:p]-Cx[(q+1:p),1:K],Wy,Wy)+
  #   NCutY3V1(Cx[(q+p+1:r),1:K],M1[1:r]-Cx[(q+p+1:r),1:K],Wx,Wx))

  Test <- vector(mode="numeric", length=B)

  for (k in 1:B){
    ###Draw k(-) and k(+)with unequal probabilites.
    if(sampling=='equal'){
      s  <- sample.int(K, K, replace=FALSE)
      while(length(which(Cx[ ,s[1]]==1))==1){
        s  <- sample.int(K, K, replace=FALSE)
      }
      ax <- which(Cx[ ,s[1]]==1)
    }else if(sampling=='size'){
      Nx <- apply(Cx, 2, sum)
      N  <- sum(Nx)
      P  <- Nx/N
      s  <- sample.int(K, K, replace=FALSE, prob=P)
      while(length(which(Cx[ ,s[1]]==1))==1){
        s  <- sample.int(K, K, replace=FALSE, prob=P)
      }
      ax <- which(Cx[ ,s[1]]==1)
    }

    sx           <- sample(ax,1)
    C2x[sx,s[1]] <- 0
    C2x[sx,s[K]] <- 1

    #Now Step 3 in the algorithm
    # J2=NCutY3V1(C2x[,1:K],M1-C2x[,1:K],Izyx2,Izyx2)+
    #   gamma*(NCutY3V1(C2x[1:q,1:K],M1[1:q]-C2x[1:q,1:K],Wz,Wz)+
    #   NCutY3V1(C2x[(q+1:p),1:K],M1[1:p]-C2x[(q+1:p),1:K],Wy,Wy)+
    #   NCutY3V1(C2x[(q+p+1:r),1:K],M1[1:r]-C2x[(q+p+1:r),1:K],Wx,Wx))

    J2 <- NCut(C2x[,1:K], Izyx2) + gamma*(NCut(C2x[1:q,1:K], Wz)+
          NCut(C2x[(q+1:p),1:K], Wy) + NCut(C2x[(q+p+1:r),1:K], Wx))

    if (J2>J){
      des=rbinom(1,1,exp(-L*log(k+1)*(J2-J)))
      if (des==1){
        Cx=C2x#Set-up the new clusters
        J=J2
        Nx=apply(Cx[,1:K],2,sum)
      }else{
        C2x=Cx
      }
    } else{
      Cx=C2x
      J=J2
      Nx=apply(Cx[,1:K],2,sum)
    }
    Test[k]=J
  }
  Res      <- list()
  Res[[1]] <- Test
  Res[[2]] <- Cx
  return(Res)
}

#' Cluster the Columns of X into K Clusters by Giving a Weighted Cluster Membership while shrinking
#' Weights Towards Each Other.
#'
#' This function will output K channels of variables.
#' @param X is a n x p matrix of p variables and n observations.
#' @param K is the number of clusters.
#' @param B is the number of iterations in the simulated annealing algorithm.
#' @param L is the temperature coefficient in the simulated annealing algorithm.
#' @param scale equals TRUE if data X is to be scaled with mean 0 and variance 1.
#' @param lambda the tuning parameter of the penalty. Larger values shrink the weighted
#' cluster membership closer together (default = 1).
#' @param nstarts the number of starting values also corresponding how many times simulated
#' annealing is run. Larger values provide better results but takes longer.
#' @param start if it equals 'default' then the starting value for all weights is 1/K. If
#' 'random' then weights are sampled from a uniform distribution and then scaled to sum 1
#' per variable.
#' @param dist specifies the distance metric used for constructing the similarity matrix.
#' Options are 'gaussian', 'correlation' and 'euclidean' (default = 'gaussian').
#' @param epsilon values in the similarity matrix less than epsilon are set to 0 (default = 0).
#' @param sigma is the bandwidth parameter when the dist metric chosen is 'gaussian' (default = 0.1).
#' @param beta when dist='correlation', beta is the exponent applied to each entry of the
#' similarity matrix.
#' @details
#' The algorithm minimizes a modified version of NCut through simulated annealing.
#' The clusters correspond to partitions that minimize this objective function.
#' @references Sebastian J. Teran Hidalgo, Mengyun Wu and Shuangge Ma.
#' Penalized and weighted clustering of gene expression data using PWNCut. (Submitted.)
#' @examples
#' # This sets up the initial parameters for the simulation.
#' n <- 100 # Sample size
#' p <- 100 # Number of columns of Y.
#' K <- 3
#'
#' C0            <- matrix(0,p,K)
#' C0[1:25,1]    <- matrix(1,25,1)
#' C0[26:75,1:3] <- matrix(1/3,50,3)
#' C0[76:100,3]  <- matrix(1,25,1)
#'
#' A0 <- C0[ ,1]%*%t(C0[ ,1]) + C0[ ,2]%*%t(C0[ ,2]) +
#'       C0[ ,3]%*%t(C0[ ,3])
#' A0 <- A0 - diag(diag(A0)) + diag(p)
#'
#' Z1 <- rnorm(n,0,2)
#' Z2 <- rnorm(n,0,2)
#' Z3 <- rnorm(n,0,2)
#'
#' Y <- matrix(0,n,p)
#' Y[ ,1:25]   <-  matrix(rnorm(n*25, 0, 2), n, 25) + matrix(Z1, n, 25, byrow=FALSE)
#' Y[ ,26:75]  <-  matrix(rnorm(n*50, 0, 2), n, 50) + matrix(Z1, n, 50, byrow=FALSE) +
#'                 matrix(Z2, n, 50, byrow=FALSE) + matrix(Z3, n, 50, byrow=FALSE)
#' Y[ ,76:100] <-  matrix(rnorm(n*25, 0, 2), n, 25) + matrix(Z3, n, 25, byrow=FALSE)
#'
#' trial <- pwncut(Y,
#'                 K       = 3,
#'                 B       = 10000,
#'                 L       = 1000,
#'                 lambda  = 1.5,
#'                 start   = 'default',
#'                 scale   = TRUE,
#'                 nstarts = 1,
#'                 epsilon = 0,
#'                 dist    = 'correlation',
#'                 sigma   = 10)
#'
#' A1 <- trial[[2]][ ,1]%*%t(trial[[2]][ ,1]) +
#'       trial[[2]][ ,2]%*%t(trial[[2]][ ,2]) +
#'       trial[[2]][ ,3]%*%t(trial[[2]][ ,3])
#'
#' A1 <- A1 - diag(diag(A1)) + diag(p)
#'
#' plot(trial[[1]], type='l')
#' errorL <- sum(abs(A0-A1))/p^2
#' errorL
#' @export
pwncut <- function(X,
                   K       = 2,
                   B       = 3000,
                   L       = 1000,
                   scale   = TRUE,
                   lambda  = 1,
                   epsilon = 0,
                   nstarts = 3,
                   start   = 'default',
                   dist    = 'gaussian',
                   sigma   = 0.1,
                   beta    = 1){
  #Beginning of the function
  if (scale==T){
    X=apply(X,2,function(e) {return(e/stats::var(e)^0.5)})
    p=dim(X)[2]
    if (dist=='gaussian'){
      Wx=exp((-1)*as.matrix(dist(t(X),diag=T,upper=T))/(sigma^2))
      Wx[which(Wx<epsilon)]=0
    }else if(dist=='euclidean'){
      Wx=as.matrix(dist(t(X),diag=T,upper=T))+diag(p)
      Wx=Wx^(-1)
      Wx[which(Wx<epsilon)]=0
    }else if(dist=='correlation'){
      Wx<-abs(stats::cor(X))^beta
      Wx[which(Wx<epsilon)]=0
    }else{
      print('Distance Error')
    }

  }else{
    p=dim(X)[2]
    if (dist=='gaussian'){
      Wx=exp((-1)*as.matrix(dist(t(X),diag=T,upper=T))/(sigma^2))
      Wx[which(Wx<epsilon)]=0
    }else if(dist=='euclidean'){
      Wx=as.matrix(dist(t(X),diag=T,upper=T))+diag(p)
      Wx=Wx^(-1)
      Wx[which(Wx<epsilon)]=0
    }else if(dist=='correlation'){
      Wx<-abs(stats::cor(X))^beta
      Wx[which(Wx<epsilon)]=0
    }else{
      print('Distance Error')
    }
  }

  #WHAT IS THIS BELOW????
  #Probs <- matrix(apply(Wx,1,sum),p,p)
  #Wx <- Wx/Probs

  #L<-diag(apply(Wx,1,sum))-Wx
  D<-matrix(apply(Wx,1,sum),p,1)
  Rfinal <- list()
  for (j in 1:nstarts){
    #This creates a random starting point in the split in the algorithm for K clusters
    if(is.matrix(start)==T){
      Cx=matrix(0,p,K)
      Cx<-start
    } else if(start=='default'){
      Cx <- matrix(1/K,p,K)
    } else if (start=='random'){
      Cx=matrix(runif(p*K),p,K)
      Sums=apply(Cx,1,sum)
      Cx=Cx/Sums
    } else if(start=='kmeans'){
      results<-kmeans(t(X),centers=K)
      Cx=matrix(0,p,K)
      for (i in 1:p){
        Cx[i,results[[1]][i]]=1
      }
    } else{
      print('This is not a start option')
    }

    #Now, calculate the number of indices in each group.
    Nx   <- apply(Cx,2,sum)
    M1   <- matrix(1,p,K)
    J    <- WNCut(Cx,M1-Cx,Wx)+lambda*Ranking(Cx)/(p*K)
    Test <- vector(mode="numeric", length=B)
      #These matrices will keep track of the elements of the clusters while
      #doing simulated annealing.
      C2x=matrix(0,p,K)
      C2x=Cx
      for (k in 1:B){
        ###New sketchs of the beginning of the algorithm###
        #1.- Draw the one of the genes randomly {i}#
        s<-sample.int(p,1)
        #2.- Rank the weights for the {i}-th gene and keep the rankings {j}#
        s.order<-order(Cx[s,])
        #3.- Draw a cluster at random#
        cluster<-sample.int(K,1)
        #4.- Decide whether make that sparse with respect to w_{i,(j-1)} or w_{i,(j+1)}#
        if(s.order[cluster]==1){#if the weight is the smallest one
          w.replace<-which(s.order==2)
          min.dist<-abs(C2x[s,s.order[cluster]]-C2x[s,w.replace])
        }else if (s.order[cluster]==K){#if the weight is the largest one
          w.replace<-which(s.order==K-1)
          min.dist<-abs(C2x[s,s.order[cluster]]-C2x[s,w.replace])
        }else{
          w1.replace<-which(s.order==(s.order[cluster]-1))
          w2.replace<-which(s.order==(s.order[cluster]+1))
          w.replace<-c(w1.replace,w2.replace)
          diffs<-c(abs(C2x[s,s.order[cluster]]-C2x[s,w1.replace]),
                   abs(C2x[s,s.order[cluster]]-C2x[s,w2.replace]))
          min.dist<-min(diffs)
          a1<-which(diffs==min.dist)
          w.replace<-w.replace[a1[1]]#what weight are we going to replace
        }

        p_minus=runif(1,min=0,max=C2x[s,cluster])
        C2x[s,cluster]=C2x[s,cluster]-p_minus#This element will give somethin between 0 and its value
        C2x[s,w.replace]=C2x[s,w.replace]+p_minus#This element will get something between 0 and the value of the other element

        #Now Step 3 in the algorithm
        J2 <- WNCut(C2x,M1-C2x,Wx) + lambda*Ranking(C2x)/(K*p)

        if (J2>J){
          des=rbinom(1,1,exp(-L*log(k+1)*(J2-J)))
          if (des==1){
            Cx=C2x#Set-up the new clusters
            J=J2
            Nx=apply(Cx,2,sum)
          }else{
            C2x=Cx
          }
        } else{
          Cx=C2x
          J=J2
          Nx=apply(Cx,2,sum)
        }
        Test[k]=J
      }
      if(j==1){
        Rfinal[[1]] <- Test
        Rfinal[[2]] <- Cx
      }else{
        if(Rfinal[[1]][B]>Test[B]){
          Rfinal[[1]] <- Test
          Rfinal[[2]] <- Cx
        }
      }
  }
  return(Rfinal)
}

#' The MLBNCut Clusters the Columns and the Rows Simultaneously of Data from 3 Different Sources.
#'
#' It clusters the columns of Z,Y and X into K clusters and the samples into R clusters by representing
#' each data type as one network layer.
#' It represents the Z layer depending on Y, and the Y layer depending on X.
#'
#' This function will output K clusters of columns of Z, Y and X and R clusters of the samples.
#' @param Z is a n x q matrix of q variables and n observations.
#' @param Y is a n x p matrix of p variables and n observations.
#' @param X is a n x r matrix of r variables and n observations.
#' @param K is the number of column clusters.
#' @param R is the number of row clusters.
#' @param B is the number of iterations.
#' @param N is the number of samples per iterations.
#' @param q0 is the quantiles in the cross entropy method.
#' @param sigmas is the tuning parameter of the Gaussian kernel of the samples.
#' @param sigmac is the tuning parameter of the Gaussian kernel of the variables.
#' @param dist is the type of distance measure use in the similarity matrix.
#' Options are 'gaussian' and 'correlation', with 'gaussian' being the default.
#' @param scale equals TRUE if data Y is to be scaled with mean 0 and variance 1.
#' @return A list with the final value of the objective function and
#' the clusters.
#' @details
#' The algorithm minimizes the NCut through the cross entropy method.
#' The clusters correspond to partitions that minimize this objective function.
#' @references Sebastian J. Teran Hidalgo and Shuangge Ma.
#' Multilayer Biclustering of Omics Data using MLBNCut. (Work in progress.)
#' @examples
#' #This sets up the initial parameters for the simulation.
#' library(NCutYX)
#' library(MASS)
#' library(fields)
#'
#' n   <- 50
#' p   <- 50
#' h   <- 0.15
#' rho <- 0.15
#' mu  <- 1
#'
#' W0 <- matrix(1,p,p)
#' W0[1:(p/5),1:(p/5)] <- 0
#' W0[(p/5+1):(3*p/5),(p/5+1):(3*p/5)] <- 0
#' W0[(3*p/5+1):(4*p/5),(3*p/5+1):(4*p/5)] <- 0
#' W0[(4*p/5+1):p,(4*p/5+1):p]=0
#' W0=cbind(W0,W0,W0)
#' W0=rbind(W0,W0,W0)
#'
#' W1 <- matrix(1,n,n)
#' W1[1:(n/2),1:(n/2)] <- 0
#' W1[(n/2+1):n,(n/2+1):n] <- 0
#'
#' X <- matrix(0,n,p)
#' Y <- matrix(0,n,p)
#' Z <- matrix(0,n,p)
#'
#' Sigma=matrix(0,p,p)
#' Sigma[1:(p/5),1:(p/5)] <- rho
#' Sigma[(p/5+1):(3*p/5),(p/5+1):(3*p/5)] <- rho
#' Sigma[(3*p/5+1):(4*p/5),(3*p/5+1):(4*p/5)] <- rho
#' Sigma[(4*p/5+1):p,(4*p/5+1):p] <- rho
#' Sigma <- Sigma - diag(diag(Sigma))
#' Sigma <- Sigma + diag(p)
#'
#' X[1:(n/2),]   <- mvrnorm(n/2,rep(mu,p),Sigma)
#' X[(n/2+1):n,] <- mvrnorm(n/2,rep(-mu,p),Sigma)
#'
#' B11 <- matrix(0,p,p)
#' B12 <- matrix(0,p,p)
#' B21 <- matrix(0,p,p)
#' B22 <- matrix(0,p,p)
#'
#' B11[1:(p/5),1:(p/5)]                     <- runif((p/5)^2,h/2,h)*rbinom((p/5)^2,1,0.5)
#' B11[(p/5+1):(3*p/5),(p/5+1):(3*p/5)]     <- runif((2*p/5)^2,h/2,h)*rbinom((2*p/5)^2,1,0.5)
#' B11[(3*p/5+1):(4*p/5),(3*p/5+1):(4*p/5)] <- runif((p/5)^2,h/2,h)*rbinom((p/5)^2,1,0.5)
#' B11[(4*p/5+1):p,(4*p/5+1):p]             <- runif((1*p/5)^2,h/2,h)*rbinom((1*p/5)^2,1,0.5)
#'
#' B12[1:(p/5),1:(p/5)]                     <- runif((p/5)^2,-h,-h/2)*rbinom((p/5)^2,1,0.5)
#' B12[(p/5+1):(3*p/5),(p/5+1):(3*p/5)]     <- runif((2*p/5)^2,-h,-h/2)*rbinom((2*p/5)^2,1,0.5)
#' B12[(3*p/5+1):(4*p/5),(3*p/5+1):(4*p/5)] <- runif((p/5)^2,-h,-h/2)*rbinom((p/5)^2,1,0.5)
#' B12[(4*p/5+1):p,(4*p/5+1):p]             <- runif((1*p/5)^2,-h,-h/2)*rbinom((1*p/5)^2,1,0.5)
#'
#' B21[1:(p/5),1:(p/5)]                     <- runif((p/5)^2,h/2,h)*rbinom((p/5)^2,1,0.5)
#' B21[(p/5+1):(3*p/5),(p/5+1):(3*p/5)]     <- runif((2*p/5)^2,h/2,h)*rbinom((2*p/5)^2,1,0.5)
#' B21[(3*p/5+1):(4*p/5),(3*p/5+1):(4*p/5)] <- runif((p/5)^2,h/2,h)*rbinom((p/5)^2,1,0.5)
#' B21[(4*p/5+1):p,(4*p/5+1):p]             <- runif((1*p/5)^2,h/2,h)*rbinom((1*p/5)^2,1,0.5)
#'
#' B22[1:(p/5),1:(p/5)]                     <- runif((p/5)^2,-h,-h/2)*rbinom((p/5)^2,1,0.5)
#' B22[(p/5+1):(3*p/5),(p/5+1):(3*p/5)]     <- runif((2*p/5)^2,-h,-h/2)*rbinom((2*p/5)^2,1,0.5)
#' B22[(3*p/5+1):(4*p/5),(3*p/5+1):(4*p/5)] <- runif((p/5)^2,-h,-h/2)*rbinom((p/5)^2,1,0.5)
#' B22[(4*p/5+1):p,(4*p/5+1):p]             <- runif((1*p/5)^2,-h,-h/2)*rbinom((1*p/5)^2,1,0.5)
#'
#' Y[1:(n/2),]   <- X[1:(n/2),]%*%B11+matrix(rnorm((n/2)*p,0,0.25),n/2,p)
#' Y[(n/2+1):n,] <- X[(n/2+1):n,]%*%B12+matrix(rnorm((n/2)*p,0,0.25),n/2,p)
#'
#' Z[1:(n/2),]   <- Y[1:(n/2),]%*%B21+matrix(rnorm((n/2)*p,0,0.25),n/2,p)
#' Z[(n/2+1):n,] <- Y[(n/2+1):n,]%*%B22+matrix(rnorm((n/2)*p,0,0.25),n/2,p)
#'
#' trial <- mlbncut(Z,
#'                  Y,
#'                  X,
#'                  K=4,
#'                  R=2,
#'                  B=10,
#'                  N=50,
#'                  dist='correlation',
#'                  q0=0.15,
#'                  scale=TRUE,
#'                  sigmas=0.05,
#'                  sigmac=1)
#'
#' plot(trial[[1]],type='l')
#' image.plot(trial[[2]])
#' image.plot(trial[[3]])
#'
#' errorK <- sum((trial[[3]][,1]%*%t(trial[[3]][,1]) +
#'                  trial[[3]][,2]%*%t(trial[[3]][,2]) +
#'                  trial[[3]][,3]%*%t(trial[[3]][,3]) +
#'                  trial[[3]][,4]%*%t(trial[[3]][,4]))*W0)/(3*p)^2 +
#'             sum((trial[[2]][,1]%*%t(trial[[2]][,1]) +
#'                  trial[[2]][,2]%*%t(trial[[2]][,2]))*W1)/(n)^2
#' errorK
#' @export
mlbncut <- function(Z,
                Y,
                X,
                K      = 2,
                R      = 2,
                B      = 30,
                N      = 500,
                q0     = 0.25,
                scale  = TRUE,
                dist   = 'gaussian',
                sigmas = 1,
                sigmac = 1){
  Res       <- list()
  quantiles <- vector(mode="numeric", length=B)
  #Beginning of the function
  if (scale==T){
    Z <- scale(Z)
    Y <- scale(Y)
    X <- scale(X)
  }
  ZYX  <- cbind(Z,Y,X)
  dimz <- dim(Z)[2]
  dimy <- dim(Y)[2]
  dimx <- dim(X)[2]
  n    <- dim(X)[1]
  m    <- dimz + dimy + dimx
  #vector with the probabilities for clustering samples and columns
  #Initialize step in the algorithm
  Ps     <- matrix(1/R,n,R)
  Pc     <- matrix(1/K,m,K)
  Clusts <- vector('list',N)
  Clustc <- vector('list',N)
  Wr     <- vector('list',R)
  Wk     <- vector('list',K)
  loss   <- vector(mode="numeric", length=N)
  #Start of the Cross entropy optimization
  #For j in {B}
  for (j in 1:B){
    print(paste('jth Loop is ', j))
    #For k in {N}
    for (k in 1:N){
      #Sample the partitions V_b^{(t)} and S_b^{(t)} from P_v^{(t-1)} and P_s^{(t-1)}
      Clustc[[k]] <- RandomMatrix(m,K,Pc)
      Clusts[[k]] <- RandomMatrix(n,R,Ps)
      loss[k]     <- 0
      if (dist=='correlation'){
        for (r in 1:R){
          c1      <- which(Clusts[[k]][,r]==1)
          Wr[[r]] <- w.cor(Z[c1, ],Y[c1, ],X[c1, ])
        }
      }

      if (dist=='gaussian'){
        for (r in 1:R){
          c1      <- which(Clusts[[k]][,r]==1)
          Wr[[r]] <- w.gauss(Z[c1, ],Y[c1, ],X[c1, ],sigma=sigmac)
        }
      }
      for (i in 1:K){
        cz      <- which(Clustc[[k]][1:dimz,i]==1)
        cy      <- which(Clustc[[k]][(dimz+1):(dimz+dimy),i]==1)
        cx      <- which(Clustc[[k]][(dimz+dimy+1):m,i]==1)
        A1      <- cbind(Z[ ,cz],Y[ ,cy],X[ ,cx])
        Wk[[i]] <- exp((-sigmas)*as.matrix(stats::dist(A1,
                                                       method = 'euclidean',
                                                       diag   = T,
                                                       upper  = T)))
      }
      for (r in 1:R){
          loss[k] <- loss[k]+NCut(Clustc[[k]],Wr[[r]])
      }
      for (i in 1:K){
          loss[k] <- loss[k]+NCut(Clusts[[k]],Wk[[i]])
      }
    }
    #Calculate \hat{q}
    cutoff       <- as.numeric(quantile(loss,q0,na.rm=T,type=4))
    quantiles[j] <- cutoff
    #Calculate P_v^{(t)} and P_s^{(t)}
    s1   <- which(loss<=cutoff)
    sumc <- Reduce('+',Clustc[s1])
    Pc   <- sumc/length(s1)
    sums <- Reduce('+',Clusts[s1])
    Ps   <- sums/length(s1)
  }#End of cross entropy optimization

  Res[[1]] <- quantiles
  Res[[2]] <- Ps
  Res[[3]] <- Pc
  return(Res)
}

#' Cluster the Rows of X into K Clusters Using the AWNCut Method.
#'
#' Builds similarity matrices for the rows of X and the rows of an assisted dataset Z.
#' Clusters them into K groups while conducting feature selection based on the AWNCut method.
#'
#' @param X is an n x p1 matrix of n observations and p1 variables.
#' @param Z is an n x p2 matrix of n observations and p2 variables. Z is the assistant dataset.
#' @param K is the number of clusters.
#' @param lambda is a vector of tuning parameter lambda in the objective function.
#' @param Tau is a vector of tuning parameters tau to be used in the objective function.
#' @param B is the number of iterations in the simulated annealing algorithm.
#' @param L is the temperature coefficient in the simulated annealing algorithm.
#' @return  A list with the following components:
#' \describe{
#' \item{lambda}{the value of tuning parameter lambda for the result}
#' \item{tau}{the value of tuning parameter tau for the result}
#' \item{Cs}{a matrix of the clustering result}
#' \item{ws}{a vector of the feature selection result}
#' \item{OP.value}{the value of the objective function}
#' }
#' @details
#' The algorithm maximizes a sum of the weighed NCut measure for X and assisted dataset Z,
#' with the addition of a correlation measure between the two datasets. Feature selection
#' is implemented by using the average correlation of each feature as a criterion.
#' @author Ruofan Bie. Maintainer: Sebastian Jose Teran Hidalgo
#' \url{sebastianteranhidalgo@gmail.com}.
#' @references Li, Yang; Bie, Ruofan; Teran Hidalgo, Sebastian; Qin, Yinchen; Wu, Mengyun; Ma, Shuangge.
#' Assisted gene expression-based clustering with AWNCut. (Submitted.)
#' @examples
#' set.seed(123456)
#' #This sets up the initial parameters for the simulation.
#' lambda <- seq(2,6,1) #Tuning parameter lambda
#' Tau    <- seq(0.2,0.8,0.2) #Tuning parameter tau
#'
#' n=30; n1=10; n2=10; n3=n-n1-n2 #Sample size
#' p1=10; p2=10; r1=8; r2=8; #Number of variables and noises in each dataset
#'
#' K=3; #Number of clusters
#'
#' mu=1; #Mean of the marginal distribution
#' u1=0.5; #Range of enties in the coefficient matrix
#'
#' library(mvtnorm)
#' epsilon <- matrix(rnorm(n*(p1-r1),0,1), n, (p1-r1)) # Generation of random error
#'
#' Sigma1 <- matrix(rep(0.8,(p1-r1)^2),(p1-r1),(p1-r1)) # Generation of the covariance matrix
#' diag(Sigma1) <- 1
#'
#' # Generation of the original distribution of the three clusters
#' T1 <- matrix(rmvnorm(n1,mean=rep(-mu,(p1-r1)),sigma=Sigma1),n1,(p1-r1))
#' T2 <- matrix(rmvnorm(n2,mean=rep(0,(p1-r1)),sigma=Sigma1),n2,(p1-r1))
#' T3 <- matrix(rmvnorm(n3,mean=rep(mu,(p1-r1)),sigma=Sigma1),n3,(p1-r1))
#'
#' X1 <- sign(T1)*(exp(abs(T1))) #Generation of signals in X
#' X2 <- sign(T2)*(exp(abs(T2)))
#' X3 <- sign(T3)*(exp(abs(T3)))
#' ep1 <- (matrix(rnorm(n*r1,0,1),n,r1)) #Generation of noises in X
#' X <- rbind(X1,X2,X3)
#'
#' beta1 <- matrix(runif((p1-r1)*(p2-r2),-u1,u1),(p1-r1),(p2-r2)) #Generation of the coefficient matrix
#' Z     <- X%*%beta1+epsilon #Generation of signals in Z
#' ep2   <- (matrix(rnorm(n*r2,0.5,1),n,r2)) #Generation of noises in Z
#'
#' X <- cbind(X,ep1)
#' Z <- cbind(Z,ep2)
#' #our method
#' Tune1         <- awncut.selection(X, Z, K, lambda, Tau, B = 20, L = 1000)
#' awncut.result <- awncut(X, Z, 3, Tune1$lam, Tune1$tau, B = 20, L = 1000)
#' ErrorRate(awncut.result[[1]]$Cs, n1, n2)
#' @export
awncut <- function(X,
                   Z,
                   K,
                   lambda,
                   Tau,
                   B = 500,
                   L = 1000){
  X <- scale(X)
  Z <- scale(Z)
  #Generate a Cartesian product of the two tuning parameters and try all possible conbinations
  Para <- as.data.frame(cbind(rep(lambda,each=length(Tau)),rep(Tau,length(lambda))))
  out  <- list()
  for(para in 1:nrow(Para)){
    #Initialization
    lam    <- Para[para,1]
    tau    <- Para[para,2]
    p1     <- ncol(X)
    p2     <- ncol(Z)
    w1     <- rep(1/sqrt(p1), p1)
    w2     <- rep(1/sqrt(p2), p2)
    b      <- 0
    ws.old <- c(w1,w2)
    ws     <- rep(0, p1+p2)
    Cs.old <- matrix(rep(0,nrow(Z)*K),nrow(Z),K)
    for(i in 1:nrow(Z)){
      Cs.old[i,sample(K,1)] <- 1
    }

    while((b<=B)||(sum(ws-ws.old)/sum(ws.old)>=10e-4)){
      b <- b+1
      #Calculate the weight datasets
      wm1 <- AWNcut.W(X, Z, ws.old)
      WX1 <- wm1[[1]]
      WZ1 <- wm1[[2]]

      #Compute the value of the objective function using the old clustering and feature selection results
      a1           <- AWNcut.OP(X, Z, WX1, WZ1, Cs.old, tau)
      OP.value.old <- a1$TOP+lam*sum(ws.old*a1$Cor.perfeature)/(p1+p2)

      #Update the clustering and feature selection results
      Cs <- AWNcut.UpdateCs(WX1, WZ1, K, Cs.old)
      ws <- AWNcut.UpdateWs(X, Z, K, WX1, WZ1, b, Cs, ws.old, tau)

      #Calculate the weight datasets using updated weights
      wm2 <- AWNcut.W(X, Z, ws)
      WX2 <- wm2[[1]]
      WZ2 <- wm2[[2]]

      #Compute the value of the objective function using the updated clustering and feature selection results
      a2       <- AWNcut.OP(X, Z, WX2, WZ2, Cs, tau)
      OP.value <- a2$TOP + lam*sum(ws*a2$Cor.perfeature)/(p1 + p2)
      if(OP.value<=OP.value.old){
        des <- rbinom(1, 1, Prob(OP.value, OP.value.old, L, b))
        if(des==1){
          Cs.old <- Cs
          ws.old <- ws
        }else{
          Cs <- Cs.old
          ws <- ws.old
        }
      }else{
        Cs.old <- Cs
        ws.old <- ws
      }
    }
    out[[para]] <- list(lambda   = lam,
                        tau      = tau,
                        Cs       = Cs.old,
                        ws       = ws.old,
                        OP.value = OP.value)
  }
  return(out)
}

#' This Function Outputs the Selection of Tuning Parameters for the AWNCut Method.
#'
#' @param X is an n x p1 matrix of n observations and p1 variables.
#' @param Z is an n x p2 matrix of n observations and p2 variables. Z is the assistant dataset.
#' @param K is the number of clusters.
#' @param lambda is a vector of tuning parameter lambda in the objective function.
#' @param Tau is a vector of tuning parameter tau in the objective function.
#' @param B is the number of iterations in the simulated annealing algorithm.
#' @param L is the temperature coefficient in the simulated annealing algorithm.
#' #' @return  A list with the following components:
#' \describe{
#' \item{num}{is the position of the max DBI}
#' \item{Table}{is the Table of the DBI for all possible combination of the parameters}
#' \item{lam}{is the best choice of tuning parameter lambda}
#' \item{tau}{is the best choice of tuning parameter lambda}
#' \item{DBI}{is the max DBI}
#' }
#' @references Li, Yang; Bie, Ruofan; Teran Hidalgo, Sebastian; Qin, Yinchen; Wu, Mengyun; Ma, Shuangge.
#' Assisted gene expression-based clustering with AWNCut. (Submitted.)
#' @examples
#' set.seed(123456)
#' #This sets up the initial parameters for the simulation.
#' lambda <- seq(2,6,1) #Tuning parameter lambda
#' Tau    <- seq(0.2,0.8,0.2) #Tuning parameter tau
#'
#' n=30; n1=10; n2=10; n3=n-n1-n2 #Sample size
#' p1=10; p2=10; r1=8; r2=8; #Number of variables and noises in each dataset
#'
#' K=3; #Number of clusters
#'
#' mu=1; #Mean of the marginal distribution
#' u1=0.5; #Range of enties in the coefficient matrix
#'
#' library(mvtnorm)
#' epsilon <- matrix(rnorm(n*(p1-r1),0,1), n, (p1-r1)) # Generation of random error
#'
#' Sigma1 <- matrix(rep(0.8,(p1-r1)^2),(p1-r1),(p1-r1)) # Generation of the covariance matrix
#' diag(Sigma1) <- 1
#'
#' # Generation of the original distribution of the three clusters
#' T1 <- matrix(rmvnorm(n1,mean=rep(-mu,(p1-r1)),sigma=Sigma1),n1,(p1-r1))
#' T2 <- matrix(rmvnorm(n2,mean=rep(0,(p1-r1)),sigma=Sigma1),n2,(p1-r1))
#' T3 <- matrix(rmvnorm(n3,mean=rep(mu,(p1-r1)),sigma=Sigma1),n3,(p1-r1))
#'
#' X1 <- sign(T1)*(exp(abs(T1))) #Generation of signals in X
#' X2 <- sign(T2)*(exp(abs(T2)))
#' X3 <- sign(T3)*(exp(abs(T3)))
#' ep1 <- (matrix(rnorm(n*r1,0,1),n,r1)) #Generation of noises in X
#' X <- rbind(X1,X2,X3)
#'
#' beta1 <- matrix(runif((p1-r1)*(p2-r2),-u1,u1),(p1-r1),(p2-r2)) #Generation of the coefficient matrix
#' Z     <- X%*%beta1+epsilon #Generation of signals in Z
#' ep2   <- (matrix(rnorm(n*r2,0.5,1),n,r2)) #Generation of noises in Z
#'
#' X <- cbind(X,ep1)
#' Z <- cbind(Z,ep2)
#' #our method
#' Tune1         <- awncut.selection(X, Z, K, lambda, Tau, B = 20, L = 1000)
#' awncut.result <- awncut(X, Z, 3, Tune1$lam, Tune1$tau, B = 20, L = 1000)
#' ErrorRate(awncut.result[[1]]$Cs, n1, n2)
#' @export
awncut.selection <- function(X,
                             Z,
                             K,
                             lambda,
                             Tau,
                             B = 500,
                             L = 1000){
  out  <- awncut(X, Z, K, lambda, Tau, B, L=1000)
  Para <- as.data.frame(cbind(rep(lambda,each=length(Tau)),rep(Tau,length(lambda))))
  dbi  <- NULL
  for(i in 1:nrow(Para)){
    Cs  <- out[[i]]$Cs
    ws  <- out[[i]]$ws
    dbi <- c(dbi, DBI(cbind(X, Z), K, Cs, ws))
  }
  return(list(num   = which.max(dbi),
              Table = t(cbind(Para, dbi)),
              lam   = Para[which.max(dbi), 1],
              tau   = Para[which.max(dbi), 2],
              DBI   = max(dbi)))
}

#' This Function Calculates the True Error Rate of a Clustering Result,
#' Assuming that There are Three Clusters.
#'
#' @return err is the true error rate of a clustering result.
#'
#' @param X is a clustering result in matrix format.
#' @param n1 is the size of the first cluster.
#' @param n2 is the size of the second cluster.
#' @references Li, Yang; Bie, Ruofan; Teran Hidalgo, Sebastian; Qin, Yinchen; Wu, Mengyun; Ma, Shuangge.
#' Assisted gene expression-based clustering with AWNCut. (Submitted.)
#' @examples
#' set.seed(123456)
#' #This sets up the initial parameters for the simulation.
#' lambda <- seq(2,6,1) #Tuning parameter lambda
#' Tau    <- seq(0.2,0.8,0.2) #Tuning parameter tau
#'
#' n=30; n1=10; n2=10; n3=n-n1-n2 #Sample size
#' p1=10; p2=10; r1=8; r2=8; #Number of variables and noises in each dataset
#'
#' K=3; #Number of clusters
#'
#' mu=1; #Mean of the marginal distribution
#' u1=0.5; #Range of enties in the coefficient matrix
#'
#' library(mvtnorm)
#' epsilon <- matrix(rnorm(n*(p1-r1),0,1), n, (p1-r1)) # Generation of random error
#'
#' Sigma1 <- matrix(rep(0.8,(p1-r1)^2),(p1-r1),(p1-r1)) # Generation of the covariance matrix
#' diag(Sigma1) <- 1
#'
#' # Generation of the original distribution of the three clusters
#' T1 <- matrix(rmvnorm(n1,mean=rep(-mu,(p1-r1)),sigma=Sigma1),n1,(p1-r1))
#' T2 <- matrix(rmvnorm(n2,mean=rep(0,(p1-r1)),sigma=Sigma1),n2,(p1-r1))
#' T3 <- matrix(rmvnorm(n3,mean=rep(mu,(p1-r1)),sigma=Sigma1),n3,(p1-r1))
#'
#' X1 <- sign(T1)*(exp(abs(T1))) #Generation of signals in X
#' X2 <- sign(T2)*(exp(abs(T2)))
#' X3 <- sign(T3)*(exp(abs(T3)))
#' ep1 <- (matrix(rnorm(n*r1,0,1),n,r1)) #Generation of noises in X
#' X <- rbind(X1,X2,X3)
#'
#' beta1 <- matrix(runif((p1-r1)*(p2-r2),-u1,u1),(p1-r1),(p2-r2)) #Generation of the coefficient matrix
#' Z     <- X%*%beta1+epsilon #Generation of signals in Z
#' ep2   <- (matrix(rnorm(n*r2,0.5,1),n,r2)) #Generation of noises in Z
#'
#' X <- cbind(X,ep1)
#' Z <- cbind(Z,ep2)
#' #our method
#' Tune1         <- awncut.selection(X, Z, K, lambda, Tau, B = 20, L = 1000)
#' awncut.result <- awncut(X, Z, 3, Tune1$lam, Tune1$tau, B = 20, L = 1000)
#' ErrorRate(awncut.result[[1]]$Cs, n1, n2)
#' @export
ErrorRate <- function(X, n1, n2){
  n <- nrow(X)
  Error <- matrix(1,n,n)
  Error[1:n1,1:n1]<-0
  Error[(1+n1):(n1+n2),(1+n1):(n1+n2)]<-0
  Error[(1+n2+n1):n,(1+n2+n1):n] <- 0
  Denum <- sum(Error)
  err <- 0
  for(i in 1:ncol(X)){
    f <- matrix(X[,i],n,1)
    err <- err+sum(f%*%t(f)*Error)/Denum
  }
  return(err)
}
