#include <RcppEigen.h>
#include <Rcpp.h>
#include <algorithm> //std::min

using namespace Rcpp;
using namespace Eigen;

//PROTOTYPE FUNCTIONS
//This is the function prototype for the overload version of NCutY3V1
double NCutY3V1(const MatrixXd &Cys, const MatrixXd &Cy2s,
                const MatrixXd &Wys, const MatrixXd &Wxs);
//Overloaded version to get Eigen matrices
MatrixXd Penal(const MatrixXd &Cys);

//This function calculates the ANCut objective function
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double NCutY3V1(const NumericMatrix &Cys, const NumericMatrix &Cy2s,
                const NumericMatrix &Wys, const NumericMatrix &Wxs){

  Eigen::Map<Eigen::MatrixXd> Wy = as<Eigen::Map<Eigen::MatrixXd> >(Wys);
  Eigen::Map<Eigen::MatrixXd> Wx = as<Eigen::Map<Eigen::MatrixXd> >(Wxs);
  Eigen::Map<Eigen::MatrixXd> Cy = as<Eigen::Map<Eigen::MatrixXd> >(Cys);
  Eigen::Map<Eigen::MatrixXd> Cy2 = as<Eigen::Map<Eigen::MatrixXd> >(Cy2s);
  const int K = Cy.cols();//number of groups
  //These vectors wil contain the number of vertices in from X and Y
  //For the K=2 case, cut1=cut2
  VectorXd Cuty(K);//Numerator of NCut
  VectorXd Volx(K);//Denoinator of NCut
  double J;
  for(int i=0;i<K;i++){
    Cuty(i)=Cy.col(i).transpose()*Wy*Cy2.col(i);//This is new
    Volx(i)=Cy.col(i).transpose()*Wx*Cy.col(i);
  }

  VectorXd Res=Cuty.array()/Volx.array();
  J=Res.sum();

  //In the general case, there needs to be  vectors
  //Volx and Voly of dimension K
  //Vol1x=(D1*Wy*D1).sum();
  //Vol2x=(D2*Wy*D2).sum();

  return J;
}

//This function calculates a NCut for a three layered graph.
//It has as input the cluster membership, and the respective
//weigth matrices.


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double NCutLayer3V1(const NumericMatrix &Cys, const NumericMatrix &Cy2s,
                    const NumericMatrix &Wzs, const NumericMatrix &Wys,
                    const NumericMatrix &Wxs, const NumericMatrix &Wzyxs){

  Eigen::Map<Eigen::MatrixXd> Wz = as<Eigen::Map<Eigen::MatrixXd> >(Wzs);
  Eigen::Map<Eigen::MatrixXd> Wy = as<Eigen::Map<Eigen::MatrixXd> >(Wys);
  Eigen::Map<Eigen::MatrixXd> Wx = as<Eigen::Map<Eigen::MatrixXd> >(Wxs);
  Eigen::Map<Eigen::MatrixXd> Wzyx = as<Eigen::Map<Eigen::MatrixXd> >(Wzyxs);
  Eigen::Map<Eigen::MatrixXd> Cy = as<Eigen::Map<Eigen::MatrixXd> >(Cys);
  Eigen::Map<Eigen::MatrixXd> Cy2 = as<Eigen::Map<Eigen::MatrixXd> >(Cy2s);

  const int q = Wz.rows();//Number of Z's
  const int p = Wy.rows();//number of Y's
  const int r = Wx.rows();//Number of X's
  const int K = Cy.cols();//number of groups

  //These vectors wil contain the number of vertices in from X and Y
  //For the K=2 case, cut1=cut2
  VectorXd Cuty(K);//Numerator of NCut
  double J;
  double a0,a1,a2,a3;
  for(int i=0;i<K;i++){
    a0=Cy.col(i).transpose()*Wzyx*Cy2.col(i);//This is new.
    a1=std::max(double(Cy.col(i).segment(0,q).transpose()*Wz*Cy.col(i).segment(0,q)),double(1));
    a2=std::max(double(Cy.col(i).segment(q,p).transpose()*Wy*Cy.col(i).segment(q,p)),double(1));
    a3=std::max(double(Cy.col(i).segment(q+p,r).transpose()*Wx*Cy.col(i).segment(q+p,r)),double(1));
    Cuty(i)=a0/(a1+a2+a3);//Why is this likes this?
  }

  J=Cuty.sum();

  //In the general case, there needs to be  vectors
  //Volx and Voly of dimension K
  //Vol1x=(D1*Wy*D1).sum();
  //Vol2x=(D2*Wy*D2).sum();

  return J;
  //return J;
}

//This function calculates a max penalty



// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
NumericMatrix Penal(const NumericMatrix &Cys){
  //Changing it to a Eigen::MatrixXd
  Eigen::Map<Eigen::MatrixXd> Cy = as<Eigen::Map<Eigen::MatrixXd> >(Cys);
  const int K = Cy.cols();//number of groups
  const int p = Cy.rows();//Number of genes
  //Penalty matrix P of penalties
  MatrixXd Pen(p,K);
  //Finding the penalties
  for (int i=0;i<p;i++){
    Pen(i,0)=(1-Cy.block(i,1,1,K-1).maxCoeff())/(K-1);
    Pen(i,K-1)=(1-Cy.block(i,0,1,K-1).maxCoeff())/(K-1);
    if (K>2){
      for (int j=1;j<K-1;j++){
        Pen(i,j)=(1-std::max(double(Cy.block(i,0,1,j).maxCoeff()),double(Cy.block(i,j+1,1,K-1-j).maxCoeff())))/(K-1);
      }
    }
  }

  NumericMatrix Penx(wrap(Pen));
  return Penx;
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::export]]
SEXP SimuA(const NumericMatrix &C0s, const NumericMatrix &Wxs,
           const double &lambda, const int &B, const double &mu_0,
           const double &L){
  //For the random number generation
  RNGScope scope;

  //Changing from NumericMatrix to MatrixXd
  Eigen::Map<Eigen::MatrixXd> C0 = as<Eigen::Map<Eigen::MatrixXd> >(C0s);
  Eigen::Map<Eigen::MatrixXd> Wx = as<Eigen::Map<Eigen::MatrixXd> >(Wxs);
  const int K = C0.rows();//number of groups
  const int p = C0.cols();//Number of genes
  double pd=p;
  double Kd=K;
  //Now, calculate the number of indices in each group.
  //Below needs to be added
  //Nx=apply(Cx,2,sum)

  //START AGAIN HERE
  MatrixXd M1 = MatrixXd::Constant(p, K, 1.0);

  //These matrices will keep track of the elements of the clusters while
  //doing simulated annealing.
  MatrixXd C2x =C0;

  double J=NCutY3V1(C0,M1-C0,Wx,Wx)+lambda*(C0-Penal(C0)).lpNorm<1>()/(K*p^2);

  VectorXd Test(B);

  for (int k=0;k<B;k++){

    int s=ceil(R::runif(0,Kd));
    int ax=ceil(R::runif(0,pd));
    //Select a vertex from cluster s[1] with unequal probability
    //Calculating Unequal probabilites

    //Draw a coin to see whether we choose X or Y
    //ax=which(Cx[,s[1]]!=0)//which Xs belong to the cluster

    //New version
    //sx=sample(ax,1)
    double Prob=mu_0+C2x(ax-1,s-1);
    int sparse=R::rbinom(1,Prob);
    //CONTINUE HERE NEXT TIME. BELOW USE BLOCK OPERATIONS AND REDUCTIONS FOR MAX
    if (sparse==0){
      p_minus=C2x[sx,s[1]]-(1-max(C2x[sx,s[2:K]]))/(K-1)
      C2x[sx,s[1]]=(1-max(C2x[sx,s[2:K]]))/(K-1)
      C2x[sx,s[2:K]]=((K-1)/K)*C2x[sx,s[2:K]]/sum(C2x[sx,s[2:K]])
    }else{
      p_minus=runif(1,min=0,max=C2x[sx,s[1]])
      C2x[sx,s[1]]=C2x[sx,s[1]]-p_minus//This element will give somethin between 0 and its value
      C2x[sx,s[K]]=C2x[sx,s[K]]+p_minus//This element will get something between 0 and the value of the other element
    }

    //Now Step 3 in the algorithm
    J2=NCutY3V1(C2x,M1-C2x,Wx,Wx)+lambda*sum(abs(C2x-Penal(C2x)))/(K*p^2)

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
      //Res<-list()
      //Res[[1]]=Test
      //Res[[2]]=Cx
      //return(Res)
    //return Rcpp::List::create(soma,div,menor,dif);
}

//OVERLOADED VERSIONS TO BE USED WITHIN C++ AND NOT R ARE HERE:

//Overloaded version of NCutY3V1 to be used within the C++ function above
//and that has a prototype function above.
double NCutY3V1(const MatrixXd &Cy, const MatrixXd &Cy2,
                const MatrixXd &Wy, const MatrixXd &Wx){

  //Eigen::Map<Eigen::MatrixXd> Wy = as<Eigen::Map<Eigen::MatrixXd> >(Wys);
  //Eigen::Map<Eigen::MatrixXd> Wx = as<Eigen::Map<Eigen::MatrixXd> >(Wxs);
  //Eigen::Map<Eigen::MatrixXd> Cy = as<Eigen::Map<Eigen::MatrixXd> >(Cys);
  //Eigen::Map<Eigen::MatrixXd> Cy2 = as<Eigen::Map<Eigen::MatrixXd> >(Cy2s);
  const int K = Cy.cols();//number of groups
  //These vectors wil contain the number of vertices in from X and Y
  //For the K=2 case, cut1=cut2
  VectorXd Cuty(K);//Numerator of NCut
  VectorXd Volx(K);//Denoinator of NCut
  double J;
  for(int i=0;i<K;i++){
    Cuty(i)=Cy.col(i).transpose()*Wy*Cy2.col(i);//This is new
    Volx(i)=Cy.col(i).transpose()*Wx*Cy.col(i);
  }

  VectorXd Res=Cuty.array()/Volx.array();
  J=Res.sum();

  //In the general case, there needs to be  vectors
  //Volx and Voly of dimension K
  //Vol1x=(D1*Wy*D1).sum();
  //Vol2x=(D2*Wy*D2).sum();

  return J;
}

//Overloaded version of the Penal function
//the prototype function is above in the prototype section
MatrixXd Penal(const MatrixXd &Cy){
  const int K = Cy.cols();//number of groups
  const int p = Cy.rows();//Number of genes
  //Penalty matrix P of penalties
  MatrixXd Pen(p,K);
  //Finding the penalties
  for (int i=0;i<p;i++){
    Pen(i,0)=(1-Cy.block(i,1,1,K-1).maxCoeff())/(K-1);
    Pen(i,K-1)=(1-Cy.block(i,0,1,K-1).maxCoeff())/(K-1);
    if (K>2){
      for (int j=1;j<K-1;j++){
        Pen(i,j)=(1-std::max(double(Cy.block(i,0,1,j).maxCoeff()),double(Cy.block(i,j+1,1,K-1-j).maxCoeff())))/(K-1);
      }
    }
  }

  return Pen;
}

