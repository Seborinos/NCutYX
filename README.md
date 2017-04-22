# Description

The NCutYX package includes functions on clustering genomic data using graph theory. 

* The current version contains the function ANCut that clusters one type of data genomic data, say gene expressions, with the help of a second type of data, like copy number aberrations. 
* The function LNCut clusters a three-layered graph into K different channels of 3 types of of genomic data. 

To install:

* latest development version: 
    1. install and load package devtools
    1. `install_github("Seborinos/NCutYX")`

# Example

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