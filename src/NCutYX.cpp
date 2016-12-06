#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

//This function calculates the NCutYX objective function.

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double NCutY3V1(const NumericMatrix Cys,const NumericMatrix Wys,
                const NumericMatrix Wxs){

  Eigen::Map<Eigen::MatrixXd> Wy = as<Eigen::Map<Eigen::MatrixXd> >(Wys);
  Eigen::Map<Eigen::MatrixXd> Wx = as<Eigen::Map<Eigen::MatrixXd> >(Wxs);
  Eigen::Map<Eigen::MatrixXd> Cy = as<Eigen::Map<Eigen::MatrixXd> >(Cys);
  const int K = Cy.cols();//number of groups
  const int p = Cy.rows();
  //These vectors wil contain the number of vertices in from X and Y
  //For the K=2 case, cut1=cut2
  VectorXd V=VectorXd::Zero(p);//need to define p
  V=V.array()+1;//vector of only ones
  VectorXd Cuty(K);//Numerator of NCut
  VectorXd Volx(K);//Denoinator of NCut
  double J;
  for(int i=0;i<K;i++){
    Cuty(i)=(Cy.col(i).transpose()*Wy*(V-Cy.col(i)))(0);
    Volx(i)=Cy.col(i).transpose()*Wx*Cy.col(i);
  }

  J=(Cuty.array()/Volx.array()).sum();

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
double NCutLayer3V1(const NumericMatrix Cys, const NumericMatrix Wzs
                      ,const NumericMatrix Wys, const NumericMatrix Wxs
                      ,const NumericMatrix Wzyxs){

  Eigen::Map<Eigen::MatrixXd> Wz = as<Eigen::Map<Eigen::MatrixXd> >(Wzs);
  Eigen::Map<Eigen::MatrixXd> Wy = as<Eigen::Map<Eigen::MatrixXd> >(Wys);
  Eigen::Map<Eigen::MatrixXd> Wx = as<Eigen::Map<Eigen::MatrixXd> >(Wxs);
  Eigen::Map<Eigen::MatrixXd> Wzyx = as<Eigen::Map<Eigen::MatrixXd> >(Wzyxs);
  Eigen::Map<Eigen::MatrixXd> Cy = as<Eigen::Map<Eigen::MatrixXd> >(Cys);

  const int q = Wz.rows();//Number of Z's
  const int p = Wy.rows();//number of Y's
  const int r = Wx.rows();//Number of X's
  const int K = Cy.cols();//number of groups

  //These vectors wil contain the number of vertices in from X and Y
  //For the K=2 case, cut1=cut2
  VectorXd V=VectorXd::Zero(p);//need to define p
  V=V.array()+1;//vector of only ones
  VectorXd Cuty(K);//Numerator of NCut
  VectorXd Volx(K);//Denoinator of NCut
  double J;
  double a1,a2,a3;
  for(int i=0;i<K;i++){
    Cuty(i)=(Cy.col(i).transpose()*Wzyx*(V-Cy.col(i)))(0);
    a1=Cy.col(i).segment(0,q).transpose()*Wz*Cy.col(i).segment(0,q);//not sure this will work
    a2=Cy.col(i).segment(q,p).transpose()*Wy*Cy.col(i).segment(q,p);//not sure about this
    a3=Cy.col(i).segment(q+p,r).transpose()*Wx*Cy.col(i).segment(q+p,r);
    Volx(3,i)=pow(a1*a2*a3,1/3);
  }

  J=(Cuty.array()/Volx.array()).sum();

  //In the general case, there needs to be  vectors
  //Volx and Voly of dimension K
  //Vol1x=(D1*Wy*D1).sum();
  //Vol2x=(D2*Wy*D2).sum();

  return J;
}
