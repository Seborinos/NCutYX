#include <RcppEigen.h>
#include <algorithm>  // std::sort std::min

using namespace Rcpp;
using namespace Eigen;

//This function calculates the NCutYX objective function.

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

//This function calculates the WNCut

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double WNCut(const NumericMatrix &Cys,
              const NumericMatrix &Cy2s,
              const NumericMatrix &Wys){

  Eigen::Map<Eigen::MatrixXd> Cy = as<Eigen::Map<Eigen::MatrixXd> >(Cys);
  Eigen::Map<Eigen::MatrixXd> Cy2 = as<Eigen::Map<Eigen::MatrixXd> >(Cy2s);
  Eigen::Map<Eigen::MatrixXd> Wy = as<Eigen::Map<Eigen::MatrixXd> >(Wys);

  const int K = Cy.cols();//number of groups
  //These vectors wil contain the number of vertices in from X and Y
  //For the K=2 case, cut1=cut2
  VectorXd Cuty(K);//Numerator of NCut
  VectorXd Volx(K);//Denoinator of NCut

  double J;
  for(int i=0;i<K;i++){
    //Cuty(i)=(Cy.col(i).transpose()*Wy*Cy2.col(i)).sum();
    Cuty(i)=(Cy.col(i).transpose()*Wy*Cy2.col(i)).sum();
    Volx(i)=(Cy.col(i).transpose()*Wy*Cy.col(i));
    //Volx(i)=(Cy.col(i).transpose()*Wy).sum();
  }

  VectorXd Res=Cuty.array()/Volx.array();
  J=Res.sum();

  //In the general case, there needs to be  vectors
  //Volx and Voly of dimension K
  //Vol1x=(D1*Wy*D1).sum();
  //Vol2x=(D2*Wy*D2).sum();

  return J;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double WNCut2(const NumericMatrix &Cys,
              const NumericMatrix &Cy2s,
              const NumericMatrix &Dys,
              const NumericMatrix &Wys){

  Eigen::Map<Eigen::MatrixXd> Cy = as<Eigen::Map<Eigen::MatrixXd> >(Cys);
  Eigen::Map<Eigen::MatrixXd> Cy2 = as<Eigen::Map<Eigen::MatrixXd> >(Cy2s);
  Eigen::Map<Eigen::MatrixXd> Dy = as<Eigen::Map<Eigen::MatrixXd> >(Dys);
  Eigen::Map<Eigen::MatrixXd> Wy = as<Eigen::Map<Eigen::MatrixXd> >(Wys);

  const int K = Cy.cols();//number of groups
  //These vectors wil contain the number of vertices in from X and Y
  //For the K=2 case, cut1=cut2
  VectorXd Cuty(K);//Numerator of NCut
  VectorXd Volx(K);//Denoinator of NCut

  double J;
  for(int i=0;i<K;i++){
    //Cuty(i)=(Cy.col(i).transpose()*Wy*Cy2.col(i)).sum();
    Cuty(i)=(Cy.col(i).transpose()*Dy*Cy2.col(i)).sum();
    Volx(i)=(Cy.col(i).transpose()*Wy*Cy.col(i));
    //Volx(i)=(Cy.col(i).transpose()*Wy).sum();
  }

  VectorXd Res=Cuty.array()/Volx.array();
  J=Res.sum();

  //In the general case, there needs to be  vectors
  //Volx and Voly of dimension K
  //Vol1x=(D1*Wy*D1).sum();
  //Vol2x=(D2*Wy*D2).sum();

  return J;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double WNCut3(const NumericMatrix &Cys,
              const NumericMatrix &Cy2s,
              const NumericMatrix &Dys,
              const NumericMatrix &Wys){

  Eigen::Map<Eigen::MatrixXd> Cy = as<Eigen::Map<Eigen::MatrixXd> >(Cys);
  Eigen::Map<Eigen::MatrixXd> Cy2 = as<Eigen::Map<Eigen::MatrixXd> >(Cy2s);
  Eigen::Map<Eigen::MatrixXd> Dy = as<Eigen::Map<Eigen::MatrixXd> >(Dys);
  Eigen::Map<Eigen::MatrixXd> Wy = as<Eigen::Map<Eigen::MatrixXd> >(Wys);

  const int K = Cy.cols();//number of groups
  //These vectors wil contain the number of vertices in from X and Y
  //For the K=2 case, cut1=cut2
  VectorXd Cuty(K);//Numerator of NCut
  VectorXd Volx(K);//Denoinator of NCut

  double J;
  for(int i=0;i<K;i++){
    //Cuty(i)=(Cy.col(i).transpose()*Wy*Cy2.col(i)).sum();
    Cuty(i)=(Cy.col(i).transpose()*Dy*Cy2.col(i)).sum();
    Volx(i)=(Cy.col(i).transpose()*Wy).sum();
    //Volx(i)=(Cy.col(i).transpose()*Wy).sum();
  }

  VectorXd Res=Cuty.array()/Volx.array();
  J=Res.sum();

  //In the general case, there needs to be  vectors
  //Volx and Voly of dimension K
  //Vol1x=(D1*Wy*D1).sum();
  //Vol2x=(D2*Wy*D2).sum();

  return J;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double WNCut4(const NumericMatrix &Cys,
              const NumericMatrix &Cy2s,
              const NumericMatrix &Wys){

  Eigen::Map<Eigen::MatrixXd> Cy = as<Eigen::Map<Eigen::MatrixXd> >(Cys);
  Eigen::Map<Eigen::MatrixXd> Cy2 = as<Eigen::Map<Eigen::MatrixXd> >(Cy2s);
  Eigen::Map<Eigen::MatrixXd> Wy = as<Eigen::Map<Eigen::MatrixXd> >(Wys);

  const int K = Cy.cols();//number of groups
  //These vectors wil contain the number of vertices in from X and Y
  //For the K=2 case, cut1=cut2
  VectorXd Cuty(K);//Numerator of NCut
  VectorXd Volx(K);//Denoinator of NCut

  double J;
  for(int i=0;i<K;i++){
    //Cuty(i)=(Cy.col(i).transpose()*Wy*Cy2.col(i)).sum();
    Cuty(i)=(Cy.col(i).transpose()*Wy*Cy2.col(i)).sum();
    Volx(i)=(Cy.col(i).transpose()*Wy).sum();
    //Volx(i)=(Cy.col(i).transpose()*Wy).sum();
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

//This function calculates a ranking penalty

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double Ranking(const NumericMatrix &C){
  //Getting the dimension of the matrix of clusters weights
  const int K = C.ncol();//number of clusters
  const int p = C.nrow();//Number of genes

  double sum=0;
  for (int i=0;i<p;i++){
    NumericVector W = C(i,_);
    std::sort (W.begin(), W.begin()+K);
    //sum=sum+W(0);
    for (int k=1;k<K;k++){
      sum=sum+(W(k)-W(k-1));
    }
  }

  return sum;
}

