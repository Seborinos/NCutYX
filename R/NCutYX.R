#' Cluster the columns of Y into 2 groups with the help of external features in X.
#'
#' This function will output 2 groups of variables.
#' @param Y is a n x p matrix of p variables and n observations.
#' @param X is a n x q matrix of q variables and n observations.
#' @param B is the number of iterations in the simulated annealing algorithm.
#' @param L is the temperature coefficient in the simulated annealing algorithm.
#' @param alpha is the coefficient of the elastic net penalty.
#' @param nlambdas is the number of tunning parameters in the elastic net.
#' @param ncv is the number of corss-validations in the elastic net.
#' @return A list with the objective function, the clusters and the lambda penalty.
#' @details
#' The algorithm minimizes a modified version of NCut through simulated annealing.
#' The clusers correspond to partitions that minimize this objective function.
#' The external information of X is incorporated by using ridge regression to predict Y.
#' This is a new comment.
NCutY2R1<-function(Y,X,B=3000,L=1000,alpha=0.5,nlambdas=100,ncv=5){
  #This creates the weight matrix W
  #W=abs(CorV1(n,p+q,cbind(X,Y)))
  #Wxy=W[1:p,(p+1):(p+q)]
  X=scale(X)
  Y=scale(Y)
  p=dim(Y)[2]
  Wyy=as.matrix(dist(t(Y),diag=T,upper=T))+diag(p)
  Wyy=Wyy^(-1)
  #This makes the any values divided by 0 into 1.
  #This should not be like this.
  #This is happening because the numbers become really small.
  #Wyy[which(Wyy==Inf)]=1

  #I changed the distance matrix to gaussian kernel
  #Wyy=exp((-1)*sigma*as.matrix(dist(t(Y),diag=T,upper=T)))

  cv.m1=cv.glmnet(X, Y, family=c("mgaussian"),
                  alpha=alpha,nfolds=ncv,nlambda=nlambdas,intercept=FALSE)
  m1=glmnet(X, Y, family=c("mgaussian"),
            alpha=alpha,lambda=cv.m1$lambda.min,intercept=FALSE)

  Y2=predict(m1,newx=X)
  Y2=scale(Y2[,,1])
  Wxx=as.matrix(dist(t(Y2),diag=T,upper=T))+diag(p)
  Wxx=Wxx^(-1)
  #Wxx[which(Wxx==Inf)]=1
  #I changed the distance matrix to gaussian kernel
  #Wxx=exp((-1)*sigma*as.matrix(dist(t(Y2),diag=T,upper=T)))
  #This creates a random starting point in the split in the algorithm for K=2
  Ax=rbinom(p,1,0.5)
  Bx=rep(1,p)-Ax

  #These matrices keep track of the elements of of the clusters 0,1
  Cx=matrix(0,p,3)
  Cx[,1]=t(Ax)
  Cx[,2]=t(Bx)
  Cx[,3]=seq(1,p,1)

  #Now, calculate the number of indices in each group.
  Nx=c()
  Nx[1]=sum(Cx[,1])
  Nx[2]=sum(Cx[,2])

  #These matrices will keep track of the elements of the clusters while
  #doing simulated annealing.
  C2x=matrix(0,p,3)
  C2x=Cx

  J=NCutY3V1(Cx[,1:2],matrix(1,p,2)-Cx[,1:2],Wyy,Wxx)

  Test<- vector(mode="numeric", length=B)

  for (k in 1:B){

    ###Draw k(-) and k(+)with unequal probabilites.
    #This section needs to change dramatically for
    #the general case

    N=Nx[1]+Nx[2]
    P=c(Nx[1]/N,Nx[2]/N)
    s=sample.int(2,2,replace=FALSE,prob=P)
    ###Select a vertex from cluster s[1] with unequal probability
    #Calculating Unequal probabilites
    #Draw a coin to see whether we choose X or Y

    ax=which(Cx[,s[1]]==1)#which Xs belong to the cluster

    sx=sample(ax,1)
    C2x[sx,s[1]]=0
    C2x[sx,s[2]]=1

    #Now Step 3 in the algorithm

    J2=NCutY2V2(C2x[,1:2],matrix(1,p,2)-C2x[,1:2],Wyy,Wxx)

    if (J2>J){
      #Prob[Count]=exp(-10000*log(k+1)*(J2-J))
      des=rbinom(1,1,exp(-L*log(k+1)*(J2-J)))
      if (des==1){
        Cx=C2x
        J=J2
        Nx[1]=sum(Cx[,1])
        Nx[2]=sum(Cx[,2])
      }else{
        C2x=Cx
      }
    } else{
      Cx=C2x
      J=J2
      Nx[1]=sum(Cx[,1])
      Nx[2]=sum(Cx[,2])
    }
    Test[k]=J

  }
  Res<-list()
  Res[[1]]=Test
  Res[[2]]=Cx
  Res[[3]]=cv.m1$lambda.min
  ##This is the return list.
  ##It contains, first...
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
#' n=200 #Sample size
#' B=5000 #Number of iterations in the simulated annealing algorithm.
#' L=1000 #Temperature coefficient.
#' p=500 #Number of columns of Y.
#' q=p #Number of columns of X.
#' h1=0
#' h2=0.15
#'
#' S=matrix(0.2,q,q)
#' S[1:(q/2),(q/2+1):q]=0
#' S[(q/2+1):q,1:(q/2)]=0
#' S=S-diag(diag(S))+diag(q)
#' S=nearPD(S)
#' S=S$mat
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
#' Res=NCutY2R1(Y,X,B,L,alpha=0,ncv=5)
#' Cx=Res[[2]]
#' f11=matrix(Cx[,1],p,1)
#' f12=matrix(Cx[,2],p,1)
#'
#' errorL=sum((f11%*%t(f11))*W0)/Denum+sum((f12%*%t(f12))*W0)/Denum
#' #This is the true error of the clustering solution.
#' errorL

NCutYR1<-function(Y,X,K=2,B=3000,L=1000,alpha=0.5,nlambdas=100,ncv=5,dist='gaussian'){
  #This creates the weight matrix W
  #W=abs(CorV1(n,p+q,cbind(X,Y)))
  #Wxy=W[1:p,(p+1):(p+q)]
  X=scale(X)
  Y=scale(Y)
  p=dim(Y)[2]
  if(dist=='euclidean'){
    Wyy=as.matrix(dist(t(Y),diag=T,upper=T))+diag(p)
    Wyy=Wyy^(-1)
  }
  if(dist=='gaussian'){
    Wyy=exp((-1)*as.matrix(dist(t(Y),diag=T,upper=T)))
  }


  #This should not be like this.
  #This is happening because the numbers become really small.
  #This needs to be changed.
  #Wyy[which(Wyy==Inf)]=1
  #I changed the distance matrix to gaussian kernel
  #Wyy=exp((-1)*sigma*as.matrix(dist(t(Y),diag=T,upper=T)))

  #This standardizes the Wxy as to make it a probability transition
  #matrix

  cv.m1=cv.glmnet(X, Y, family=c("mgaussian"),
                  alpha=alpha,nfolds=ncv,nlambda=nlambdas,intercept=FALSE)
  m1=glmnet(X, Y, family=c("mgaussian"),
            alpha=alpha,lambda=cv.m1$lambda.min,intercept=FALSE)

  Y2=predict(m1,newx=X)
  Y2=scale(Y2[,,1])
  if(dist=='euclidean'){
    Wxx=as.matrix(dist(t(Y2),diag=T,upper=T))+diag(p)
    Wxx=Wxx^(-1)
  }
  if(dist=='gaussian'){
    Wxx=exp((-1)*as.matrix(dist(t(Y2),diag=T,upper=T)))
  }

  #Wxx[which(Wxx==Inf)]=1
  #I changed the distance matrix to gaussian kernel
  #Wxx=exp((-1)*sigma*as.matrix(dist(t(Y2),diag=T,upper=T)))

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
    s=sample.int(K,K,replace=FALSE,prob=P)

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

NCutYLayer3R1<-function(Z,Y,X,K=2,B=3000,L=1000,alpha=0.5,ncv=3,nlambdas=100){
  #This creates the weight matrix W
  #W=abs(CorV1(n,p+q,cbind(X,Y)))
  #Wxy=W[1:p,(p+1):(p+q)]
  Z=scale(Z)
  Y=scale(Y)
  X=scale(X)
  q=dim(Z)[2]
  p=dim(Y)[2]
  r=dim(X)[2]
  #Joint distance matrix
  ZYX=cbind(Z,Y,X)
  Wzyx=as.matrix(dist(t(ZYX),diag=T,upper=T))+diag(q+p+r)
  Wzyx=Wzyx^(-1)
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
  #Distance matrix for the predicted variables
  ZYX2=cbind(Z2,Y2,X)
  Wzyx2=as.matrix(dist(t(ZYX2),diag=T,upper=T))+diag(q+p+r)
  Wzyx2=Wzyx2^(-1)
  #Z's distance matrix
  Wz=as.matrix(dist(t(Z2),diag=T,upper=T))+diag(q)
  Wz=Wz^(-1)
  #Y's distance matrix
  Wy=as.matrix(dist(t(Y2),diag=T,upper=T))+diag(p)
  Wy=Wy^(-1)
  #X's distance matrix
  Wx=as.matrix(dist(t(X),diag=T,upper=T))+diag(r)
  Wx=Wx^(-1)
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
  #J=NCutY3V1(Cx[,1:K],M1-Cx[,1:K],Wzyx2,Wzyx)
  J=NCutLayer3V1(Cx[,1:K],M1-Cx[,1:K],Wz,Wy,Wx,Wzyx)

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
    #J2=NCutY3V1(C2x[,1:K],M1-C2x[,1:K],Wzyx2,Wzyx)
    J2=NCutLayer3V1(C2x[,1:K],M1-C2x[,1:K],Wz,Wy,Wx,Wzyx)

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



