#include <RcppEigen.h>
#include <algorithm> //std::min

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
    a1=Cy.col(i).segment(0,q).transpose()*Wz*Cy.col(i).segment(0,q);//not sure this will work
    a2=Cy.col(i).segment(q,p).transpose()*Wy*Cy.col(i).segment(q,p);//not sure about this
    a3=Cy.col(i).segment(q+p,r).transpose()*Wx*Cy.col(i).segment(q+p,r);
    Cuty(i)=a0/pow(a1*a2*a3,0.33333333333333333333333);//Why is this likes this?
  }

  J=Cuty.sum();

  //In the general case, there needs to be  vectors
  //Volx and Voly of dimension K
  //Vol1x=(D1*Wy*D1).sum();
  //Vol2x=(D2*Wy*D2).sum();

  return J;
  //return J;
}

//prototype of overloaded function
double NCutLayer3V1(const Eigen::MatrixXd &Cy, const Eigen::MatrixXd &Cy2,
                    const Eigen::MatrixXd &Wz, const Eigen::MatrixXd &Wy,
                    const Eigen::MatrixXd &Wx, const Eigen::MatrixXd &Wzyx);

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
NumericMatrix OptimForward(const NumericMatrix &Wzs, const NumericMatrix &Wys,
                           const NumericMatrix &Wxs, const NumericMatrix &Wzyxs,
                           const int &m){
  //changing from R object to Eigen type
  Eigen::Map<Eigen::MatrixXd> Wz = as<Eigen::Map<Eigen::MatrixXd> >(Wzs);
  Eigen::Map<Eigen::MatrixXd> Wy = as<Eigen::Map<Eigen::MatrixXd> >(Wys);
  Eigen::Map<Eigen::MatrixXd> Wx = as<Eigen::Map<Eigen::MatrixXd> >(Wxs);
  Eigen::Map<Eigen::MatrixXd> Wzyx = as<Eigen::Map<Eigen::MatrixXd> >(Wzyxs);
  //Move as sequence and calculate the best sequence of graph partitions going forward.
  //Here I need to create a starting point where a first vertex is chosen.
  MatrixXd Cx=MatrixXd::Zero(m,2);
  MatrixXd Cx2=MatrixXd::Zero(m,2);
  Cx2=Cx2.array()+1;//This should work for sure
  Cx.col(1)=Cx.col(1).array()+1;//Not sure if this makes sense
  MatrixXd Res(m,2);
  for (int j=0;j<m;j++){
    int count=0;
    for (int i=0;i<m;i++){//This iterates through
      if(Cx(i,0)==0){
        Cx(i,0)=1;
        Res(i,0)=NCutLayer3V1(Cx,Cx2-Cx,Wz,Wy,Wx,Wzyx);//This calculates all of the cuts
        Cx(i,0)=0;
      }
      if(Res(i,0)<Res(count,0)){
        count=i;
        Res(i,1)=i;
      }
    }
    Cx(count,0)=1;
  }
  Rcpp::NumericMatrix Resx(wrap(Res));
  return wrap(Resx);
}

//Overloading NCutLayer3V1
double NCutLayer3V1(const Eigen::MatrixXd &Cy, const Eigen::MatrixXd &Cy2,
                    const Eigen::MatrixXd &Wz, const Eigen::MatrixXd &Wy,
                    const Eigen::MatrixXd &Wx, const Eigen::MatrixXd &Wzyx){

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
    a0=(Cy.col(i).transpose()*Wzyx*Cy2.col(i))(0,0);
    a1=std::max(double(Cy.col(i).segment(0,q).transpose()*Wz*Cy.col(i).segment(0,q)),double(1));
    a2=std::max(double(Cy.col(i).segment(q,p).transpose()*Wy*Cy.col(i).segment(q,p)),double(1));
    a3=std::max(double(Cy.col(i).segment(q+p,r).transpose()*Wx*Cy.col(i).segment(q+p,r)),double(1));
    Cuty(i)=a0/pow(a1*a2*a3,0.33333333333333333333333);//Why is this likes this?
  }

  J=Cuty.sum();

  return J;
  //return J;
}
