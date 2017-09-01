#' Cluster the columns of Y into K groups using the NCut graph measure.
#'
#' This function will output K clusters of the columns of Y.
#' @param Y is a n x p matrix of p variables and n observations. The p columns of
#' Y will be clustered into K groups using NCut.
#' @param K is the number of clusters.
#' @param B is the number of iterations.
#' @param N is the number of samples per iterations.
#' @param scale equals TRUE if data Y is to be scaled with mean 0 and variance 1.
#' @return A list with the final value of the objective function and
#' the clusters.
#' @details
#' The algorithm minimizes the NCut through the cross entropy method.
#' The clusters correspond to partitions that minimize this objective function.
#' @examples
#' #This sets up the initial parameters for the simulation.
#' library(MASS)
#' n=200 #Sample size
#' B=30 #Number of iterations in the simulated annealing algorithm.
#' p=500 #Number of columns of Y.
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
#' #Our method
#' Res=NCut(Y,B=30,N=500,K=2,dist='gaussian',sigma=1)
#' Cx=Res[[2]]
#' f11=matrix(Cx[,1],p,1)
#' f12=matrix(Cx[,2],p,1)
#'
#' errorL=sum((f11%*%t(f11))*W0)/Denum+sum((f12%*%t(f12))*W0)/Denum
#' #This is the true error of the clustering solution.
#' errorL

ncut<-function(Y,
               K=2,
               B=30,
               N=500,
               dist='correlation',
               scale=T,
               q=0.1,
               sigma=1){
  #This creates the weight matrix W
  Res <- list()
  quantiles <- vector(mode="numeric", length=B)
  if(scale==T){
    Y=scale(Y)
  }

  p=dim(Y)[2]
  if(dist=='euclidean'){
    Wyy=as.matrix(dist(t(Y),diag=T,upper=T))+diag(p)
    Wyy=Wyy^(-1)
  }else if(dist=='gaussian'){
    Wyy=exp((-1)*as.matrix(dist(t(Y),diag=T,upper=T))/sigma)
  }else if(dist=='correlation'){
    Wyy<-abs(cor(Y))
  }else{
    print('Distance Error')
  }

  #vector with probabilities of mus being 0 or not
  Ps <- matrix(1/K,p,K)
  for (j in 1:B){
    print(paste('jth Loop is ', j))
    Clusters <- vector('list',N)
    loss <- vector(mode="numeric", length=N)
    #M1 <- matrix(1L,p,K)
    for (k in 1:N){
       Clusters[[k]]=RandomMatrix(p,K,Ps)
       loss[k] <- NCut(Clusters[[k]],Wyy)#Can we make this run faster?
    }

    cutoff <- quantile(loss,q)
    s1 <- which(loss<=cutoff)
    quantiles[j] <- cutoff
    sums <- Reduce('+',Clusters[s1])
    Ps <- sums/length(s1)

  }
  Res[[1]] <- quantiles
  Res[[2]] <- Ps
  return(Res)
}

#' Cluster the columns of Y into K groups with the help of external features in X.
#'
#' This function will output K clusters of  the columns of Y, using the help of
#' X.
#' @param Y is a n x p matrix of p variables and n observations. The columns of
#' Y will be clustered into K groups.
#' @param X is a n x q matrix of q variables and n observations.
#' @param B is the number of iterations in the simulated annealing algorithm.
#' @param L is the temperature coefficient in the simulated annealing algorithm.
#' @param alpha is the coefficient of the elastic net penalty.
#' @param nlambdas is the number of tunning parameters in the elastic net.
#' @param ncv is the number of cross-validations in the elastic net.
#' @return A list with the final value of the objective function,
#' the clusters and the lambda penalty chosen through cross-validation.
#' @details
#' The algorithm minimizes a modified version of NCut through simulated annealing.
#' The clusters correspond to partitions that minimize this objective function.
#' The external information of X is incorporated by using ridge regression to predict Y.
#' @examples
#' #This sets up the initial parameters for the simulation.
#' library(MASS)#for mvrnorm
#' library(fields)#for image.plot
#' n=200 #Sample size
#' B=5000 #Number of iterations in the simulated annealing algorithm.
#' L=10000 #Temperature coefficient.
#' p=200 #Number of columns of Y.
#' q=p #Number of columns of X.
#' h1=0.05
#' h2=0.15
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
#' Y=Z+matrix(rnorm(n*p,0,2),n,p)
#' #Our method
#' Res=ancut(Y=Y,X=X,B=B,L=L,alpha=0,ncv=5)
#' Cx=Res[[2]]
#' f11=matrix(Cx[,1],p,1)
#' f12=matrix(Cx[,2],p,1)
#'
#' errorL=sum((f11%*%t(f11))*W0)/Denum+sum((f12%*%t(f12))*W0)/Denum
#' #This is the true error of the clustering solution.
#' errorL
#' #Below is a plot of the simulated annealing path.
#' plot(Res[[1]],type='l')
#' #Cluster found by ANCut
#' image.plot(Cx)

ancut<-function(Y,
                X,
                K=2,
                B=3000,
                L=1000,
                alpha=0.5,
                nlambdas=100,
                sampling='equal',
                ncv=5,
                dist='correlation',
                sigma=1){
  #This creates the weight matrix W
  X=scale(X)
  Y=scale(Y)
  p=dim(Y)[2]

  p=dim(Y)[2]
  if(dist=='euclidean'){
    Wyy=as.matrix(dist(t(Y),diag=T,upper=T))+diag(p)
    Wyy=Wyy^(-1)
  }else if(dist=='gaussian'){
    Wyy=exp((-1)*as.matrix(dist(t(Y),diag=T,upper=T))/sigma)
  }else if(dist=='correlation'){
    Wyy<-abs(cor(Y))
  }else{
    print('Distance Error')
  }

  #modelling the relationship between Y and X
  cv.m1=cv.glmnet(X, Y, family=c("mgaussian"),
                  alpha=alpha,nfolds=ncv,nlambda=nlambdas,intercept=FALSE)
  m1=glmnet(X, Y, family=c("mgaussian"),
            alpha=alpha,lambda=cv.m1$lambda.min,intercept=FALSE)

  Y2=predict(m1,newx=X)
  Y2=scale(Y2[,,1])
  if(dist=='euclidean'){
    Wxx=as.matrix(dist(t(Y2),diag=T,upper=T))+diag(p)
    Wxx=Wxx^(-1)
  }else if(dist=='gaussian'){
    Wxx=exp((-1)*as.matrix(dist(t(Y2),diag=T,upper=T))/sigma)
  }else if(dist=='correlation'){
    Wxx<-abs(cor(Y2))
  }else{
    print('Distance Error')
  }

  #This creates a random starting point in the split in the algorithm for K clusters
  Cx=matrix(0,p,K)

  for (i in 1:p){
    Cx[i,sample(K,1)]=1
  }

  #Now, calculate the number of indices in each group.
  Nx=apply(Cx[,1:K],2,sum)

  #These matrices will keep track of the elements of the clusters while
  #doing simulated annealing.
  C2x=matrix(0,p,K+1)
  C2x=Cx

  J=NCutY3V1(Cx[,1:K],matrix(1,p,K)-Cx[,1:K],Wyy,Wxx)

  Test<- vector(mode="numeric", length=B)

  for (k in 1:B){
    ###Draw k(-) and k(+)with unequal probabilites.
    #This section needs to change dramatically for
    #the general case
    N=sum(Nx)
    P=Nx/N

    if(sampling=='equal'){
      s=sample.int(K,K,replace=FALSE)
    }else if(sampling=='size'){
      s=sample.int(K,K,replace=FALSE,prob=P)
    }

    ###Select a vertex from cluster s[1] with unequal probability
    #Calculating Unequal probabilites
    #Draw a coin to see whether we choose X or Y
    ax=which(Cx[,s[1]]==1)#which Xs belong to the cluster

    sx=sample(ax,1)
    C2x[sx,s[1]]=0
    C2x[sx,s[K]]=1

    #Now Step 3 in the algorithm
    J2=NCutY3V1(C2x[,1:K],matrix(1,p,K)-C2x[,1:K],Wyy,Wxx)

    if (J2>J){
      #Prob[Count]=exp(-10000*log(k+1)*(J2-J))
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
  Res<-list()
  Res[[1]]=Test
  Res[[2]]=Cx
  Res[[3]]=cv.m1$lambda.min
  return(Res)
}

#' Cluster the columns of Z,Y and X into K channels.
#'
#' This function will output K channels of variables.
#' @param Z is a n x p1 matrix of p1 variables and n observations.
#' @param Y is a n x p2 matrix of p2 variables and n observations.
#' @param X is a n x p3 matrix of p3 variables and n observations.
#' @param B is the number of iterations in the simulated annealing algorithm.
#' @param L is the temperature coefficient in the simulated annealing algorithm.
#' @details
#' The algorithm minimizes a modified version of NCut through simulated annealing.
#' The clusers correspond to partitions that minimize this objective function.
#' The external information of X is incorporated by using ridge regression to predict Y.

muncut<-function(Z,
                 Y,
                 X,
                 K=2,
                 B=3000,
                 L=1000,
                 alpha=0.5,
                 ncv=3,
                 nlambdas=100,
                 scale=F,
                 model=F,
                 gamma=0.5,
                 dist='gaussian',
                 sigma=1){
  #Beginning of the function
  if (model==T){
    if (scale==T){
      Z=scale(Z)
      Y=scale(Y)
      X=scale(X)
      q=dim(Z)[2]
      p=dim(Y)[2]
      r=dim(X)[2]
      #Elastic net to predict Y with X
      cv.m1=cv.glmnet(X, Y, family=c("mgaussian"),
                      alpha=alpha,nfolds=ncv,nlambda=nlambdas,intercept=FALSE)
      m1=glmnet(X, Y, family=c("mgaussian"),
                alpha=alpha,lambda=cv.m1$lambda.min,intercept=FALSE)
      Y2=predict(m1,newx=X)
      Y2=scale(Y2[,,1])
      #Elastic net to predict Z with Y
      cv.m2=cv.glmnet(Y, Z, family=c("mgaussian"),
                      alpha=alpha,nfolds=ncv,nlambda=nlambdas,intercept=FALSE)
      m2=glmnet(Y, Z, family=c("mgaussian"),
                alpha=alpha,lambda=cv.m1$lambda.min,intercept=FALSE)
      Z2=predict(m2,newx=Y)
      Z2=scale(Z2[,,1])
      if (dist=='euclidean'){
        #Distance matrix for the predicted variables
        ZYX2=cbind(Z2,Y2,X)
        Wzyx2=as.matrix(dist(t(ZYX2),diag=T,upper=T))+diag(q+p+r)
        Wzyx2=Wzyx2^(-1)
        #Z's distance matrix
        Wz=as.matrix(dist(t(Z),diag=T,upper=T))+diag(q)
        Wz=Wz^(-1)
        #Y's distance matrix
        Wy=as.matrix(dist(t(Y),diag=T,upper=T))+diag(p)
        Wy=Wy^(-1)
        #X's distance matrix
        Wx=as.matrix(dist(t(X),diag=T,upper=T))+diag(r)
        Wx=Wx^(-1)
        #Matrix without diagonal entries
        Izyx2=Wzyx2
        Izyx2[1:q,1:q]=0
        Izyx2[(1:p+q),(1:p+q)]=0
        Izyx2[(1:r+p+q),(1:r+p+q)]=0
      }else if(dist=='gaussian'){
        #Distance matrix for the predicted variables
        ZYX2=cbind(Z2,Y2,X)
        Wzyx2=exp((-1)*as.matrix(dist(t(ZYX2),diag=T,upper=T))/sigma)
        #Z's distance matrix
        Wz=exp((-1)*as.matrix(dist(t(Z),diag=T,upper=T))/sigma)
        #Y's distance matrix
        Wy=exp((-1)*as.matrix(dist(t(Y),diag=T,upper=T))/sigma)
        #X's distance matrix
        Wx=exp((-1)*as.matrix(dist(t(X),diag=T,upper=T))/sigma)
        #Matrix without diagonal entries
        Izyx2=Wzyx2
        Izyx2[1:q,1:q]=0
        Izyx2[(1:p+q),(1:p+q)]=0
        Izyx2[(1:r+p+q),(1:r+p+q)]=0
      } else{
        print('Distance Error')
      }


    }else{
      q=dim(Z)[2]
      p=dim(Y)[2]
      r=dim(X)[2]
      #Elastic net to predict Y with X
      cv.m1=cv.glmnet(X, Y, family=c("mgaussian"),
                      alpha=alpha,nfolds=ncv,nlambda=nlambdas,intercept=FALSE)
      m1=glmnet(X, Y, family=c("mgaussian"),
                alpha=alpha,lambda=cv.m1$lambda.min,intercept=FALSE)
      Y2=predict(m1,newx=X)
      Y2=Y2[,,1]
      #Elastic net to predict Z with Y
      cv.m2=cv.glmnet(Y, Z, family=c("mgaussian"),
                      alpha=alpha,nfolds=ncv,nlambda=nlambdas,intercept=FALSE)
      m2=glmnet(Y, Z, family=c("mgaussian"),
                alpha=alpha,lambda=cv.m2$lambda.min,intercept=FALSE)
      Z2=predict(m2,newx=Y)
      Z2=Z2[,,1]
      if (dist=='euclidean'){
        #Distance matrix for the predicted variables
        ZYX2=cbind(Z2,Y2,X)
        Wzyx2=as.matrix(dist(t(ZYX2),diag=T,upper=T))+diag(q+p+r)
        Wzyx2=Wzyx2^(-1)
        #Z's distance matrix
        Wz=as.matrix(dist(t(Z),diag=T,upper=T))+diag(q)
        Wz=Wz^(-1)
        #Y's distance matrix
        Wy=as.matrix(dist(t(Y),diag=T,upper=T))+diag(p)
        Wy=Wy^(-1)
        #X's distance matrix
        Wx=as.matrix(dist(t(X),diag=T,upper=T))+diag(r)
        Wx=Wx^(-1)
        #Matrix without diagonal entries
        Izyx2=Wzyx2
        Izyx2[1:q,1:q]=0
        Izyx2[(1:p+q),(1:p+q)]=0
        Izyx2[(1:r+p+q),(1:r+p+q)]=0
      }else if(dist=='gaussian'){
        #Distance matrix for the predicted variables
        ZYX2=cbind(Z2,Y2,X)
        Wzyx2=exp((-1)*as.matrix(dist(t(ZYX2),diag=T,upper=T))/sigma)
        #Z's distance matrix
        Wz=exp((-1)*as.matrix(dist(t(Z),diag=T,upper=T))/sigma)
        #Y's distance matrix
        Wy=exp((-1)*as.matrix(dist(t(Y),diag=T,upper=T))/sigma)
        #X's distance matrix
        Wx=exp((-1)*as.matrix(dist(t(X),diag=T,upper=T))/sigma)
        #Matrix without diagonal entries
        Izyx2=Wzyx2
        Izyx2[1:q,1:q]=0
        Izyx2[(1:p+q),(1:p+q)]=0
        Izyx2[(1:r+p+q),(1:r+p+q)]=0
      } else{
        print('Distance Error')
      }

    }
  }else{
    if (scale==T){
      Z=scale(Z)
      Y=scale(Y)
      X=scale(X)
      q=dim(Z)[2]
      p=dim(Y)[2]
      r=dim(X)[2]
      if (dist=='euclidean'){
        #Distance matrix for the predicted variables
        ZYX2=cbind(Z,Y,X)
        Wzyx2=as.matrix(dist(t(ZYX2),diag=T,upper=T))+diag(q+p+r)
        Wzyx2=Wzyx2^(-1)
        #Z's distance matrix
        Wz=as.matrix(dist(t(Z),diag=T,upper=T))+diag(q)
        Wz=Wz^(-1)
        #Y's distance matrix
        Wy=as.matrix(dist(t(Y),diag=T,upper=T))+diag(p)
        Wy=Wy^(-1)
        #X's distance matrix
        Wx=as.matrix(dist(t(X),diag=T,upper=T))+diag(r)
        Wx=Wx^(-1)
        #Matrix without diagonal entries
        Izyx2=Wzyx2
        Izyx2[1:q,1:q]=0
        Izyx2[(1:p+q),(1:p+q)]=0
        Izyx2[(1:r+p+q),(1:r+p+q)]=0
      }else if(dist=='gaussian'){
        #Distance matrix for the predicted variables
        ZYX2=cbind(Z,Y,X)
        Wzyx2=exp((-1)*as.matrix(dist(t(ZYX2),diag=T,upper=T))/sigma)
        #Z's distance matrix
        Wz=exp((-1)*as.matrix(dist(t(Z),diag=T,upper=T))/sigma)
        #Y's distance matrix
        Wy=exp((-1)*as.matrix(dist(t(Y),diag=T,upper=T))/sigma)
        #X's distance matrix
        Wx=exp((-1)*as.matrix(dist(t(X),diag=T,upper=T))/sigma)
        #Matrix without diagonal entries
        Izyx2=Wzyx2
        Izyx2[1:q,1:q]=0
        Izyx2[(1:p+q),(1:p+q)]=0
        Izyx2[(1:r+p+q),(1:r+p+q)]=0
      } else{
        print('Distance Error')
      }

    }else{
      q=dim(Z)[2]
      p=dim(Y)[2]
      r=dim(X)[2]
      if (dist=='euclidean'){
        #Distance matrix for the predicted variables
        ZYX2=cbind(Z,Y,X)
        Wzyx2=as.matrix(dist(t(ZYX2),diag=T,upper=T))+diag(q+p+r)
        Wzyx2=Wzyx2^(-1)
        #Z's distance matrix
        Wz=as.matrix(dist(t(Z),diag=T,upper=T))+diag(q)
        Wz=Wz^(-1)
        #Y's distance matrix
        Wy=as.matrix(dist(t(Y),diag=T,upper=T))+diag(p)
        Wy=Wy^(-1)
        #X's distance matrix
        Wx=as.matrix(dist(t(X),diag=T,upper=T))+diag(r)
        Wx=Wx^(-1)
        #Matrix without diagonal entries
        Izyx2=Wzyx2
        Izyx2[1:q,1:q]=0
        Izyx2[(1:p+q),(1:p+q)]=0
        Izyx2[(1:r+p+q),(1:r+p+q)]=0
      }else if(dist=='gaussian'){
        #Distance matrix for the predicted variables
        ZYX2=cbind(Z,Y,X)
        Wzyx2=exp((-1)*as.matrix(dist(t(ZYX2),diag=T,upper=T))/sigma)
        #Z's distance matrix
        Wz=exp((-1)*as.matrix(dist(t(Z),diag=T,upper=T))/sigma)
        #Y's distance matrix
        Wy=exp((-1)*as.matrix(dist(t(Y),diag=T,upper=T))/sigma)
        #X's distance matrix
        Wx=exp((-1)*as.matrix(dist(t(X),diag=T,upper=T))/sigma)
        #Matrix without diagonal entries
        Izyx2=Wzyx2
        Izyx2[1:q,1:q]=0
        Izyx2[(1:p+q),(1:p+q)]=0
        Izyx2[(1:r+p+q),(1:r+p+q)]=0
      } else{
        print('Distance Error')
      }

    }
  }


  #This creates a random starting point in the split in the algorithm for K clusters
  Cx=matrix(0,q+p+r,K)
  #This is a matrix of only ones
  M1<-matrix(1,q+p+r,K)
  #Below: force the result to have one member per type of data and cluster.
  #The code below makes sure each cluster gets at least Min elements
  #per data type
  Min=2
  Check=matrix(0,3,K)
  while(sum((Check<=Min))>0){
    for (i in 1:(q+p+r)){
      Cx[i,sample(K,1)]=1
    }
    Check[1,]=apply(Cx[1:q,1:K],2,sum)
    Check[2,]=apply(Cx[(q+1):(q+p),1:K],2,sum)
    Check[3,]=apply(Cx[(p+q+1):(q+p+r),1:K],2,sum)
  }

  #Now, calculate the number of indices in each group.
  Nx=apply(Cx[,1:K],2,sum)
  #Nx=Check

  #These matrices will keep track of the elements of the clusters while
  #doing simulated annealing.
  C2x=matrix(0,p+q+r,K)
  C2x=Cx

  J=NCutY3V1(Cx[,1:K],M1-Cx[,1:K],Izyx2,Izyx2)+
    gamma*(NCutY3V1(Cx[1:q,1:K],M1[1:q]-Cx[1:q,1:K],Wz,Wz)+
    NCutY3V1(Cx[(q+1:p),1:K],M1[1:p]-Cx[(q+1:p),1:K],Wy,Wy)+
    NCutY3V1(Cx[(q+p+1:r),1:K],M1[1:r]-Cx[(q+p+1:r),1:K],Wx,Wx))

  Test<- vector(mode="numeric", length=B)

  for (k in 1:B){
    ###Draw k(-) and k(+)with unequal probabilites.
    #This section needs to change dramatically for
    #the general case

    N=sum(Nx)
    P=Nx/N
    s=sample.int(K,K,replace=FALSE)#No probability right now

    ###Select a vertex from cluster s[1] with unequal probability
    #Calculating Unequal probabilites
    #Draw a coin to see whether we choose X or Y
    ax=which(Cx[,s[1]]==1)#which Xs belong to the cluster

    sx=sample(ax,1)
    C2x[sx,s[1]]=0
    C2x[sx,s[K]]=1

    #Now Step 3 in the algorithm
    J2=NCutY3V1(C2x[,1:K],M1-C2x[,1:K],Izyx2,Izyx2)+
      gamma*(NCutY3V1(C2x[1:q,1:K],M1[1:q]-C2x[1:q,1:K],Wz,Wz)+
      NCutY3V1(C2x[(q+1:p),1:K],M1[1:p]-C2x[(q+1:p),1:K],Wy,Wy)+
      NCutY3V1(C2x[(q+p+1:r),1:K],M1[1:r]-C2x[(q+p+1:r),1:K],Wx,Wx))

    if (J2>J){
      #Prob[Count]=exp(-10000*log(k+1)*(J2-J))
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
  Res<-list()
  Res[[1]]=Test
  Res[[2]]=Cx
  return(Res)
}

#' Cluster the columns of X into K clusters by giving a weight for each cluster while remaining sparse.
#'
#' This function will output K channels of variables.
#' @param X is a n x p matrix of p variables and n observations.
#' @param B is the number of iterations in the simulated annealing algorithm.
#' @param L is the temperature coefficient in the simulated annealing algorithm.
#' @details
#' The algorithm minimizes a modified version of NCut through simulated annealing.
#' The clusers correspond to partitions that minimize this objective function.

spawn<-function(X,
                K=2,
                B=3000,
                L=1000,
                N=100,
                scale=T,
                lambda=1,
                epsilon=0,
                beta=1,
                nstarts=10,
                start='default',
                dist='gaussian',
                sigma=1){
  #Beginning of the function
  if (scale==T){
    X=apply(X,2,function(e) {return(e/var(e)^0.5)})
    p=dim(X)[2]
    if (dist=='gaussian'){
      Wx=exp((-1)*as.matrix(dist(t(X),diag=T,upper=T))/(sigma^2))
      Wx[which(Wx<epsilon)]=0
    }else if(dist=='euclidean'){
      Wx=as.matrix(dist(t(X),diag=T,upper=T))+diag(p)
      Wx=Wx^(-1)
      Wx[which(Wx<epsilon)]=0
    }else if(dist=='correlation'){
      Wx<-abs(cor(X))^beta
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
      Wx<-abs(cor(X))^beta
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
    Nx=apply(Cx,2,sum)
    M1=matrix(1,p,K)
    #J=WNCut3(Cx,M1-Cx,Wx)+lambda*Ranking(Cx)/(K*p)
    J=WNCut(Cx,M1-Cx,Wx)+lambda*Ranking(Cx)/(p*K)
    #J=WNCut(Cx,M1-Cx,Wx)+lambda*Ranking6(Cx,alpha)/(p*K)
    #J=WNCut(Cx,M1-Cx,Wx)+lambda*sum(abs(Cx-1/K))/(p*K)
    Test<- vector(mode="numeric", length=B)
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
        J2=WNCut(C2x,M1-C2x,Wx)+lambda*Ranking(C2x)/(K*p)
        #J2=WNCut(C2x,M1-C2x,Wx)+lambda*Ranking6(C2x,alpha)/(K*p)
        #J2=WNCut(C2x,M1-C2x,Wx)+lambda*sum(abs(C2x-1/K))/(p*K)
        #J2=WNCut3(C2x,M1-C2x,Wx)+lambda*Ranking(C2x)/(K*p)

        if (J2>J){
          #Prob[Count]=exp(-10000*log(k+1)*(J2-J))
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

#' Cluster the columns of Y into K groups using the NCut graph measure.
#'
#' This function will output K clusters of the columns of Y.
#' @param Y is a n x p matrix of p variables and n observations. The p columns of
#' Y will be clustered into K groups using NCut.
#' @param K is the number of clusters.
#' @param B is the number of iterations.
#' @param N is the number of samples per iterations.
#' @param scale equals TRUE if data Y is to be scaled with mean 0 and variance 1.
#' @return A list with the final value of the objective function and
#' the clusters.
#' @details
#' The algorithm minimizes the NCut through the cross entropy method.
#' The clusters correspond to partitions that minimize this objective function.
#' @examples
#' #This sets up the initial parameters for the simulation.
#' library(MASS)
#' n=200 #Sample size
#' B=30 #Number of iterations in the simulated annealing algorithm.
#' p=500 #Number of columns of Y.
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
#' #Our method
#' Res=NCut(Y,B=30,N=500,K=2,dist='gaussian',sigma=1)
#' Cx=Res[[2]]
#' f11=matrix(Cx[,1],p,1)
#' f12=matrix(Cx[,2],p,1)
#'
#' errorL=sum((f11%*%t(f11))*W0)/Denum+sum((f12%*%t(f12))*W0)/Denum
#' #This is the true error of the clustering solution.
#' errorL

spawn2<-function(X,
               K=2,
               B=30,
               N=500,
               dist='correlation',
               scale=T,
               q=0.1,
               beta=1,
               lambda=1,
               epsilon=0,
               sigma=1){
  #This creates the weight matrix W
  Res <- list()
  quantiles <- vector(mode="numeric", length=B)
  #Beginning of the function
  if (scale==T){
    X=apply(X,2,function(e) {return(e/var(e)^0.5)})
    p=dim(X)[2]
    if (dist=='gaussian'){
      Wx=exp((-1)*as.matrix(dist(t(X),diag=T,upper=T))/(sigma^2))
      Wx[which(Wx<epsilon)]=0
    }else if(dist=='euclidean'){
      Wx=as.matrix(dist(t(X),diag=T,upper=T))+diag(p)
      Wx=Wx^(-1)
      Wx[which(Wx<epsilon)]=0
    }else if(dist=='correlation'){
      Wx<-abs(cor(X))^beta
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
      Wx<-abs(cor(X))^beta
      Wx[which(Wx<epsilon)]=0
    }else{
      print('Distance Error')
    }
  }

  M1=matrix(1,p,K)
  #Pmin contains the minumum value of the uniform distribution for sample for the weights
  #Pmax contains the minumum value of the uniform distribution for sample for the weights
  Pmin <- matrix(0,p,K)
  Pmax <- matrix(1/K,p,K)
  for (j in 1:B){
    Clusters <- vector('list',N)
    loss <- vector(mode="numeric", length=N)

    for (k in 1:N){
      Clusters[[k]]=RandomUnifMatrix(p,K,Pmin,Pmax)
      Probs <- matrix(apply(Clusters[[k]],1,sum),p,K)
      Clusters[[k]] <- Clusters[[k]]/Probs
      loss[k] <- WNCut(Clusters[[k]],M1-Clusters[[k]],Wx)+lambda*Ranking(Clusters[[k]])/(p*K)
    }

    cutoff <- quantile(loss,q)
    s1 <- which(loss<=cutoff)
    quantiles[j] <- cutoff
    Pmin <- Reduce('matrixMIN',Clusters[s1])
    Pmax <- Reduce('matrixMAX',Clusters[s1])

  }
  Res[[1]] <- quantiles
  P <- (Pmin+Pmax)/2
  Probs <- matrix(apply(P,1,sum),p,K)
  Res[[2]] <- P/Probs
  return(Res)
}

#' Cluster the columns of Y into K groups using the NCut graph measure.
#'
#' This function will output K clusters of the columns of Y.
#' @param Y is a n x p matrix of p variables and n observations. The p columns of
#' Y will be clustered into K groups using NCut.
#' @param K is the number of clusters.
#' @param B is the number of iterations.
#' @param N is the number of samples per iterations.
#' @param scale equals TRUE if data Y is to be scaled with mean 0 and variance 1.
#' @return A list with the final value of the objective function and
#' the clusters.
#' @details
#' The algorithm minimizes the NCut through the cross entropy method.
#' The clusters correspond to partitions that minimize this objective function.
#' @examples
#' #This sets up the initial parameters for the simulation.
#' library(MASS)
#' n=200 #Sample size
#' B=30 #Number of iterations in the simulated annealing algorithm.
#' p=500 #Number of columns of Y.
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
#' #Our method
#' Res=NCut(Y,B=30,N=500,K=2,dist='gaussian',sigma=1)
#' Cx=Res[[2]]
#' f11=matrix(Cx[,1],p,1)
#' f12=matrix(Cx[,2],p,1)
#'
#' errorL=sum((f11%*%t(f11))*W0)/Denum+sum((f12%*%t(f12))*W0)/Denum
#' #This is the true error of the clustering solution.
#' errorL

# bilayer<-function(Z,Y,X
#                   K=2,
#                   R=2,
#                   B=30,
#                   N=500,
#                   dist='correlation',
#                   scale=T,
#                   q=0.1,
#                   sigma=1){
#   #This creates the weight matrix W
#   Res <- list()
#   quantiles <- vector(mode="numeric", length=B)
#   #Beginning of the function
#   if (scale==T){
#     Z=scale(Z)
#     Y=scale(Y)
#     X=scale(X)
#   }
#   q=dim(Z)[2]
#   p=dim(Y)[2]
#   r=dim(X)[2]
#   m=q+p+r
#   if (dist=='euclidean'){
#     #Distance matrix for the predicted variables
#     ZYX2=cbind(Z,Y,X)
#     Wzyx2=as.matrix(dist(t(ZYX2),diag=T,upper=T))+diag(q+p+r)
#     Wzyx2=Wzyx2^(-1)
#     Izyx2=Wzyx2
#     Izyx2[1:q,1:q]=0
#     Izyx2[(1:p+q),(1:p+q)]=0
#     Izyx2[(1:r+p+q),(1:r+p+q)]=0
#   }else if(dist=='gaussian'){
#     #Distance matrix for the predicted variables
#     ZYX2=cbind(Z,Y,X)
#     Wzyx2=exp((-1)*as.matrix(dist(t(ZYX2),diag=T,upper=T))/sigma)
#     #Matrix without diagonal entries
#     Izyx2=Wzyx2
#     Izyx2[1:q,1:q]=0
#     Izyx2[(1:p+q),(1:p+q)]=0
#     Izyx2[(1:r+p+q),(1:r+p+q)]=0
#   }else if(dist=='correlation'){
#     ZYX2=cbind(Z,Y,X)
#     Wzyx2 <- abs(cor(ZYX2))
#     #Matrix without diagonal entries
#     Izyx2=Wzyx2
#     Izyx2[1:q,1:q]=0
#     Izyx2[(1:p+q),(1:p+q)]=0
#     Izyx2[(1:r+p+q),(1:r+p+q)]=0
#   }else{
#     print('Distance Error')
#   }
#
#   #vector with the probabilites for clustering samples and columns
#   Ps <- matrix(1/R,n,R)
#   Pc <- matrix(1/K,m,K)
#   Dclust <- vector('list',R)
#   #Start of the Cross entropy optimization
#   for (j in 1:B){
#     print(paste('jth Loop is ', j))
#     Clusts <- vector('list',N)
#     Clustc <- vector('list',N)
#     loss <- vector(mode="numeric", length=N)
#     #1.-Cluster the columns
#     for (i in 1:R){
#       Dclust[[i]] <- w.gaussian(Z,Y,X,sigma2=1,w=Ps[,i])#maybe rescale the weigths to equal 1
#     }
#     #Sampling N new clusters and calculating the corresponding loss
#     for (k in 1:N){
#       Clustc[[k]]=RandomMatrix(m,K,Pc)
#       loss[k]<-0
#       for (i in 1:R){#for each of the sample groups(THIS IS A DOUBLE LOOP)
#         loss[k] <- loss[k]+NCut(Clustc[[k]],Dclust[[i]])
#       }
#     }
#
#     cutoff <- quantile(loss,q)
#     s1 <- which(loss<=cutoff)
#     sumc <- Reduce('+',Clustc[s1])
#     Pc <- sumc/length(s1)
#
#     #CONTINUE HERE
#     #CONTINUE HERE
#     #CONTINUE HERE
#     #2.-Cluster the samples
#     #Sampling N new clusters and calculating the corresponding loss
#     for (k in 1:N){
#       Clusts[[k]]=RandomMatrix(n,R,Ps)
#       for (i in 1:R){
#         Dclust[[i]] <- w.gaussian(Z,Y,X,sigma2=1,w=Clusts[[k]][,i])
#       }
#       loss[k]<-0
#       for (i in 1:R){#for each of the sample groups
#         loss[k] <- loss[k]+WNCut(Pc,Dclust[[i]])
#       }
#     }
#
#     cutoff <- quantile(loss,q)
#     s1 <- which(loss<=cutoff)
#     quantiles[j] <- cutoff
#     sums <- Reduce('+',Clusts[s1])
#     Ps <- sums/length(s1)
#
#   }#End of cross entropy optimization
#
#   Res[[1]] <- quantiles
#   Res[[2]] <- Ps
#   Res[[3]] <- Pc
#   return(Res)
# }

#' Cluster the columns of X into K nonexhaustive overlapping clusters.
#'
#' This function will output K weighted clusters of variables.
#' @param X is a n x p matrix of p variables and n observations.
#' @param B is the number of iterations in the simulated annealing algorithm.
#' @param L is the temperature coefficient in the simulated annealing algorithm.
#' @details
#' The algorithm minimizes a modified version of NCut through simulated annealing.
#' The clusers correspond to weighted partitions that minimize this objective function.
#' @examples
#' #This sets up the initial parameters for the simulation.
#' library(NCutYX)
#' K=4
#' n=300
#' p=120
#'
#' Z1<-rnorm(n,0,2)
#' Z2<-rnorm(n,0,2)
#' Z3<-rnorm(n,0,2)
#' Z4<-rnorm(n,0,2)
#'
#' X=matrix(0,n,p)
#' X[,1:20]=matrix(rnorm(n*20,0,0.5),n,20)+matrix(Z1,n,20,byrow=F)
#' X[,21:40]=matrix(rnorm(n*20,0,0.5),n,20)+0.5*matrix(Z1,n,20,byrow=F)+0.5*matrix(Z2,n,20,byrow=F)
#' X[,41:60]=matrix(rnorm(n*20,0,0.5),n,20)+matrix(Z2,n,20,byrow=F)
#' X[,61:80]=matrix(rnorm(n*20,0,0.5),n,20)+matrix(Z3,n,20,byrow=F)
#' X[,81:100]=matrix(rnorm(n*20,0,0.5),n,20)+0.5*matrix(Z3,n,20,byrow=F)+0.5*matrix(Z4,n,20,byrow=F)
#' X[,101:120]=matrix(rnorm(n*20,0,0.5),n,20)+matrix(Z4,n,20,byrow=F)
#'
#' C0<-matrix(0,p,K)
#' C0[1:20,1]<-matrix(1,1:20,1)
#' C0[21:40,1:2]<-matrix(0.5,1:20,1:2)
#' C0[41:60,2]<-matrix(1,1:20,1)
#' C0[61:80,3]<-matrix(1,1:20,1)
#' C0[81:100,3:4]<-matrix(0.5,1:20,1:2)
#' C0[101:120,4]<-matrix(1,1:20,1)
#'
#' A0 <- C0[,1]%*%t(C0[,1]) + C0[,2]%*%t(C0[,2]) +
#'     C0[,3]%*%t(C0[,3]) + C0[,4]%*%t(C0[,4])
#'
#' A0 <- A0-diag(diag(A0))+diag(p)
#'
#' result <- neoncut(X,K=2,
#'                     B=3000,
#'                     L=1000,
#'                     scale=F,
#'                     alpha=0.5,
#'                     beta=0.5,
#'                     dist='gaussian',
#'                     sigma=1)
#'
#'trial <- result[[2]]
#'
#' A <- trial[,1]%*%t(trial[,1]) + trial[,2]%*%t(trial[,2]) +
#'        trial[,3]%*%t(trial[,3]) + trial[,4]%*%t(trial[,4])
#'
#' A <- A-diag(diag(A))+diag(p)
#'
#' sum(abs(A0-A))/p^2

neoncut<-function(X,
                  K=2,
                  B=3000,
                  L=1000,
                  scale=F,
                  alpha=0.5,
                  beta=0.5,
                  dist='gaussian',
                  sigma=1){
  #Beginning of the function
  if (scale==T){
    X=scale(X)
    p=dim(X)[2]
    if (dist=='gaussian'){
      Wx=exp((-1)*as.matrix(dist(t(X),diag=T,upper=T))/sigma)
    }else if(dist=='euclidean'){
      Wx=as.matrix(dist(t(X),diag=T,upper=T))+diag(p)
      Wx=Wx^(-1)
    } else{
      print('Distance Error')
    }
  }else{
    p=dim(X)[2]
    if (dist=='gaussian'){
      Wx=exp((-1)*as.matrix(dist(t(X),diag=T,upper=T))/sigma)
    }else if(dist=='euclidean'){
      Wx=as.matrix(dist(t(X),diag=T,upper=T))+diag(p)
      Wx=Wx^(-1)
    } else{
      print('Distance Error')
    }
  }

  #This creates a random starting point in the split in the algorithm for K clusters
  Cx=matrix(rbinom(p*K,1,1/K),p,K)

  #Now, calculate the number of indices in each group.
  Nx=apply(Cx,2,sum)
  M1=matrix(1,p,K)
  #These matrices will keep track of the elements of the clusters while
  #doing simulated annealing.
  C2x=matrix(0,p,K)
  C2x=Cx
  Penal2<-length(which(apply(Cx,1,sum)==0))
  J=NCutY3V1(Cx,M1-Cx,Wx,Wx)+alpha*sum(diag(t(Cx)%*%Cx))+beta*Penal2

  Test<- vector(mode="numeric", length=B)

  for (k in 1:B){

    ###Draw k(-) and k(+)with unequal probabilites.
    N=sum(Nx)
    P=Nx/N
    s=sample.int(K,K,replace=FALSE,prob=P)
    sx=sample(p,1)
    if (C2x[sx,s[1]]==1){
      C2x[sx,s[1]]=0
    }else{
      C2x[sx,s[1]]=1
    }

    #Now Step 3 in the algorithm
    Penal2<-length(which(apply(Cx,1,sum)==0))
    J2=NCutY3V1(C2x,M1-C2x,Wx,Wx)+alpha*sum(diag(t(C2x)%*%C2x))+beta*Penal2

    if (J2>J){
      #Prob[Count]=exp(-10000*log(k+1)*(J2-J))
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
  Res<-list()
  Res[[1]]=Test
  Res[[2]]=Cx
  return(Res)
}

#' Cluster the columns of X into K nonexhaustive overlapping clusters.
#'
#' This function will output K weighted clusters of variables.
#' @param X is a n x p matrix of p variables and n observations.
#' @param B is the number of iterations in the cross entropy maximization algorithm.
#' @param N is the number of samples at each iteration.
#' @param q is the proportion of top results out of N used to recalculate the
#' sampling probabilities at each iteration.
#' @details
#' The algorithm minimizes a modified version of NCut through simulated annealing.
#' The clusers correspond to weighted partitions that minimize this objective function.
#' @examples
#' #This sets up the initial parameters for the simulation.
#' library(NCutYX)
#' K=4
#' n=300
#' p=120
#'
#' Z1<-rnorm(n,0,2)
#' Z2<-rnorm(n,0,2)
#' Z3<-rnorm(n,0,2)
#' Z4<-rnorm(n,0,2)
#'
#' X=matrix(0,n,p)
#' X[,1:20]=matrix(rnorm(n*20,0,0.5),n,20)+matrix(Z1,n,20,byrow=F)
#' X[,21:40]=matrix(rnorm(n*20,0,0.5),n,20)+
#' 0.5*matrix(Z1,n,20,byrow=F)+0.5*matrix(Z2,n,20,byrow=F)
#' X[,41:60]=matrix(rnorm(n*20,0,0.5),n,20)+matrix(Z2,n,20,byrow=F)
#' X[,61:80]=matrix(rnorm(n*20,0,0.5),n,20)+matrix(Z3,n,20,byrow=F)
#' X[,81:100]=matrix(rnorm(n*20,0,0.5),n,20)+
#' 0.5*matrix(Z3,n,20,byrow=F)+0.5*matrix(Z4,n,20,byrow=F)
#' X[,101:120]=matrix(rnorm(n*20,0,0.5),n,20)+matrix(Z4,n,20,byrow=F)
#'
#' C0<-matrix(0,p,K)
#' C0[1:20,1]<-matrix(1,1:20,1)
#' C0[21:40,1:2]<-matrix(0.5,1:20,1:2)
#' C0[41:60,2]<-matrix(1,1:20,1)
#' C0[61:80,3]<-matrix(1,1:20,1)
#' C0[81:100,3:4]<-matrix(0.5,1:20,1:2)
#' C0[101:120,4]<-matrix(1,1:20,1)
#'
#' A0 <- C0[,1]%*%t(C0[,1]) + C0[,2]%*%t(C0[,2]) +
#'     C0[,3]%*%t(C0[,3]) + C0[,4]%*%t(C0[,4])
#'
#' A0 <- A0-diag(diag(A0))+diag(p)
#'
#' result <- neoncut2(X,K=4,
#'                     B=30,
#'                     N=1000,
#'                     scale=F,
#'                     alpha=0.5,
#'                     beta=0.5,
#'                     dist='gaussian',
#'                     sigma=1)
#'
#'trial <- result[[2]]
#'
#' A <- trial[,1]%*%t(trial[,1]) + trial[,2]%*%t(trial[,2]) +
#'        trial[,3]%*%t(trial[,3]) + trial[,4]%*%t(trial[,4])
#'
#' A <- A-diag(diag(A))+diag(p)
#'
#' sum(abs(A0-A))/p^2

neoncut2<-function(X,
                  K=2,
                  B=30,
                  N=1000,
                  scale=F,
                  alpha=0.5,
                  beta=0.5,
                  epsilon=0,
                  q=0.10,
                  dist='correlation',
                  sigma=1){
  #This creates the weight matrix W
  Res <- list()
  quantiles <- vector(mode="numeric", length=B)
  #Beginning of the function
  if (scale==T){
    X=apply(X,2,function(e) {return(e/var(e)^0.5)})
    p=dim(X)[2]
    if (dist=='gaussian'){
      Wx=exp((-1)*as.matrix(dist(t(X),diag=T,upper=T))/(sigma^2))
      Wx[which(Wx<epsilon)]=0
    }else if(dist=='euclidean'){
      Wx=as.matrix(dist(t(X),diag=T,upper=T))+diag(p)
      Wx=Wx^(-1)
      Wx[which(Wx<epsilon)]=0
    }else if(dist=='correlation'){
      Wx<-abs(cor(X))^beta
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
      Wx<-abs(cor(X))^beta
      Wx[which(Wx<epsilon)]=0
    }else{
      print('Distance Error')
    }
  }

  M1=matrix(1,p,K)
  #vector with probabilities of mus being 0 or not
  Ps <- matrix(1/K,p,K)
  for (j in 1:B){
    print(paste('jth Loop is ', j))
    Clusters <- vector('list',N)
    loss <- vector(mode="numeric", length=N)

    for (k in 1:N){
      Clusters[[k]]=t(apply(Ps,1,neovector))#Can I write a faster version in Rcpp?
      Penal2<-length(which(apply(Clusters[[k]],1,sum)==0))
      loss[k] <- NCutY3V1(Clusters[[k]],M1-Clusters[[k]],Wx,Wx)+
                 alpha*sum(diag(t(Clusters[[k]])%*%Clusters[[k]]))+
                 beta*Penal2
    }

    cutoff <- quantile(loss,q)
    s1 <- which(loss<=cutoff)
    quantiles[j] <- cutoff
    sums <- Reduce('+',Clusters[s1])
    Ps <- sums/length(s1)

  }
  Res[[1]] <- quantiles
  Res[[2]] <- Ps
  return(Res)
}

#' Given a matrix U of K cluster membership and data X calculate the clusters' means.
#'
#' This function will output K clusters' means.
#' @param X is a n x p matrix of p variables and n observations.
#' @details
#' Given a matrix U of K cluster membership calculated the clusters' means.

meanx<-function(U,X){
  n<-dim(U)[1]
  K<-dim(U)[2]
  p<-dim(X)[2]
  M<-matrix(0,K,p)
  for (i in 1:K){
    s<-which(U[,i]==1)
    M[i,]<-apply(X[s,],2,mean)
  }
  return(M)
}

#' Given a data matrix X and matrix M of centroids calculate the distance
#' with respect to the centroids for each data matrix.
#'
#' This function will output a matrix K times p matrix of distances.
#' @param X is a n x p matrix of p variables and n observations.
#' @param M is a K x p matrix of cluster centroids.
#' @details
#' This is a function needed for the neokmeans algorithm.

distm<-function(M,X){
  K<-dim(M)[1]
  p<-dim(X)[2]
  D<-matrix(0,n,K)
  for (i in 1:K){
    Mx<-matrix(M[i,],n,p,byrow=T)
    Dx<-(X-Mx)^2
    D[,i]=appy(Dx,1,sum)
  }
  return(D)
}

#' Cluster the columns of X into K nonexhaustive overlapping clusters.
#'
#' This function will output K clusters of variables.
#' @param X is a n x p matrix of p variables and n observations.
#' @param B is the number of iterations in the simulated annealing algorithm.
#' @param L is the temperature coefficient in the simulated annealing algorithm.
#' @details
#' The algorithm minimizes a modified version of NCut through simulated annealing.
#' The clusers correspond to partitions that minimize this objective function.

neokmeans<-function(X,
                  K=2,
                  B=3000,
                  L=1000,
                  mu_0=0.01,
                  alpha=0.5,
                  beta=0.5,
                  MCMC=T){
  #getting the dimensions
  n<-dim(X)[1]
  p<-dim(X)[2]

  #This creates a random starting point in the split in the algorithm for K clusters
  Cx=matrix(rbinom(p*K,1,1/K),p,K)

  #Now, calculate the number of indices in each group.
  Nx=apply(Cx,2,sum)

  #These matrices will keep track of the elements of the clusters while
  #doing simulated annealing.
  C2x=matrix(0,p,K)
  C2x=Cx
  Penal2<-length(which(apply(Cx,1,sum)==0))
  Mx<-meanx(Cx,X)
  Dx<-distm(Mx,X)
  J=sum(Dx*Cx)+(1+alpha)*n*sum(diag(t(Cx)%*%Cx))+beta*n*Penal2

  Test<- vector(mode="numeric", length=B)

  for (k in 1:B){

    ###Draw k(-) and k(+)with unequal probabilites.
    N=sum(Nx)
    P=Nx/N
    s=sample.int(K,K,replace=FALSE,prob=P)
    sx=sample(p,1)
    if (C2x[sx,s[1]]==1){
      C2x[sx,s[1]]=0
    }else{
      C2x[sx,s[1]]=1
    }

    #Now Step 3 in the algorithm
    Penal2<-length(which(apply(C2x,1,sum)==0))
    M2x<-meanx(C2x,X)
    D2x<-distm(M2x,X)
    J2<-sum(D2x*C2x)+(1+alpha)*n*sum(diag(t(C2x)%*%C2x))+beta*n*Penal2

    if (J2>J){
      #Prob[Count]=exp(-10000*log(k+1)*(J2-J))
      des=rbinom(1,1,exp(-L*log(k+1)*(J2-J)))*(1-MCMC)+rbinom(1,1,1-min(c(J2/J,1)))*MCMC

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
  Res<-list()
  Res[[1]]=Test
  Res[[2]]=Cx
  return(Res)
}

