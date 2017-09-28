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

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double NCut(const NumericMatrix &Cys,
            const NumericMatrix &Wys){

  Eigen::Map<Eigen::MatrixXd> Wy = as<Eigen::Map<Eigen::MatrixXd> >(Wys);
  Eigen::Map<Eigen::MatrixXd> Cy = as<Eigen::Map<Eigen::MatrixXd> >(Cys);
  const int K = Cy.cols();//number of groups
  const int p = Cy.rows();
  MatrixXd Cy2 = MatrixXd::Ones(p,K);
  //These vectors wil contain the number of vertices in from X and Y
  //For the K=2 case, cut1=cut2
  VectorXd Cuty(K);//Numerator of NCut
  VectorXd Volx(K);//Denoinator of NCut
  double J;
  for(int i=0;i<K;i++){
    Cuty(i)=Cy.col(i).transpose()*Wy*(Cy2.col(i)-Cy.col(i));//This is new
    Volx(i)=Cy.col(i).transpose()*Wy*Cy.col(i);
  }

  VectorXd Res=Cuty.array()/Volx.array();
  J=Res.sum();

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
    Volx(i)=(Cy.col(i).transpose()*Wy*Cy.col(i)).sum();
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
    Volx(i)=(Cy.col(i).transpose()*Wy*Cy.col(i)).sum();
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

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double WNCut4(const NumericMatrix &Cys,
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
double WNCut5(const NumericMatrix &Cys,
              const NumericMatrix &Dys,
             const NumericMatrix &Wys){

  Eigen::Map<Eigen::MatrixXd> Cy = as<Eigen::Map<Eigen::MatrixXd> >(Cys);
  Eigen::Map<Eigen::MatrixXd> Dy = as<Eigen::Map<Eigen::MatrixXd> >(Dys);
  Eigen::Map<Eigen::MatrixXd> Wy = as<Eigen::Map<Eigen::MatrixXd> >(Wys);

  const int K = Cy.cols();//number of groups
  //These vectors wil contain the number of vertices in from X and Y
  //For the K=2 case, cut1=cut2
  VectorXd Cuty(K);//Numerator of NCut
  VectorXd Volx(K);//Denoinator of NCut

  double J;
  for(int i=0;i<K;i++){
    Cuty(i)=0;
    Volx(i)=(Cy.col(i).transpose()*Wy*Cy.col(i));
    for(int j=0;j<K;j++){
      Cuty(i)=(Cy.col(i).transpose()*Dy*Cy.col(j)).sum()+Cuty(i);
    }
    Cuty(i)=Cuty(i)-(Cy.col(i).transpose()*Dy*Cy.col(i)).sum();
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
double WNCut6(const NumericMatrix &Cys,
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
    Volx(i)=Cy.col(i).sum();
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
double WNCut7(const NumericMatrix &Cys,
              const NumericMatrix &Cy2s,
              const NumericMatrix &Dys){

  Eigen::Map<Eigen::MatrixXd> Cy = as<Eigen::Map<Eigen::MatrixXd> >(Cys);
  Eigen::Map<Eigen::MatrixXd> Cy2 = as<Eigen::Map<Eigen::MatrixXd> >(Cy2s);
  Eigen::Map<Eigen::MatrixXd> Dy = as<Eigen::Map<Eigen::MatrixXd> >(Dys);

  const int K = Cy.cols();//number of groups
  //These vectors wil contain the number of vertices in from X and Y
  //For the K=2 case, cut1=cut2
  VectorXd Cuty(K);//Numerator of NCut
  VectorXd Volx(K);//Denoinator of NCut

  double J;
  for(int i=0;i<K;i++){
    //Cuty(i)=(Cy.col(i).transpose()*Wy*Cy2.col(i)).sum();
    Cuty(i)=(Cy.col(i).transpose()*Dy*Cy2.col(i)).sum();
    Volx(i)=Cy.col(i).sum();
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
double WNCut8(const NumericMatrix &Cys,
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
    Cuty(i)=(Cy.col(i).transpose()*Wy).sum();
    Volx(i)=Cy.col(i).sum();
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
      sum=sum+(W(k)-W(k-1))*(W(k)-W(k-1));
    }
  }

  return sum;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double Ranking2(const NumericMatrix &C){
  //Getting the dimension of the matrix of clusters weights
  const int K = C.ncol();//number of clusters
  const int p = C.nrow();//Number of genes

  double sum=0;
  for (int i=0;i<p;i++){
    NumericVector W = C(i,_);
    std::sort (W.begin(), W.begin()+K);
    for (int k=0;k<K-1;k++){
      for (int j =k+1;j<K;j++){
        sum=sum+(W(j)-W(k))*(W(j)-W(k));
      }
    }
  }

  return sum;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double Ranking3(const NumericMatrix &C){
  //Getting the dimension of the matrix of clusters weights
  const int K = C.ncol();//number of clusters
  const int p = C.nrow();//Number of genes

  double sum=0;
  for (int i=0;i<p;i++){
    NumericVector W = C(i,_);
    std::sort (W.begin(), W.begin()+K);
      sum=sum+(W(K-1)-W(0));
  }
  return sum;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double Ranking4(const NumericMatrix &C){
  //Getting the dimension of the matrix of clusters weights
  const int K = C.ncol();//number of clusters
  const int p = C.nrow();//Number of genes

  double sum=0;
  for (int i=0;i<p;i++){
    NumericVector W = C(i,_);
    std::sort (W.begin(), W.begin()+K);
    sum=sum+(W(K-1)-W(0))*(W(K-1)-W(0));
    for (int k=1;k<K;k++){
      sum=sum+(W(k)-W(k-1))*(W(k)-W(k-1));
    }
  }

  return sum;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double Ranking5(const NumericMatrix &C){
  //Getting the dimension of the matrix of clusters weights
  const int K = C.ncol();//number of clusters
  const int p = C.nrow();//Number of genes

  double sum=0;
  for (int i=0;i<p;i++){
    NumericVector W = C(i,_);
    std::sort (W.begin(), W.begin()+K);
    sum=sum+(W(K-1)-W(0));
    for (int k=1;k<K;k++){
      sum=sum+(W(k)-W(k-1));
    }
  }

  return sum;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double Ranking6(const NumericMatrix &C, const double alpha){
  //Getting the dimension of the matrix of clusters weights
  const int K = C.ncol();//number of clusters
  const int p = C.nrow();//Number of genes

  double sum=0;
  for (int i=0;i<p;i++){
    NumericVector W = C(i,_);
    std::sort (W.begin(), W.begin()+K);
    //sum=sum+W(0);
    for (int k=1;k<K;k++){
      sum=sum+alpha*(W(k)-W(k-1))*(W(k)-W(k-1))+(1-alpha)*(W(k)-W(k-1));
    }
  }

  return sum;
}


// [[Rcpp::export]]
IntegerVector oneMultinomCalt(NumericVector probs) {
  int k = probs.size();
  IntegerVector ans(k);
  rmultinom(1, probs.begin(), k, ans.begin());
  return(ans);
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
IntegerMatrix RandomMatrix(const int &p,
                           const int &K,
                           const NumericMatrix &P){
  IntegerMatrix Cluster(p,K);
  //std::fill(Cluster.begin(), Cluster.end(),0);

  for (int i=0;i<p;i++){
        Cluster(i,_)=oneMultinomCalt(P(i,_));
  }

  return Cluster;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
NumericMatrix RandomUnifMatrix(const int &p,
                               const int &K,
                               const NumericMatrix &Pmin,
                               const NumericMatrix &Pmax){
  NumericMatrix Cluster(p,K);
  //std::fill(Cluster.begin(), Cluster.end(),0);

  for (int i=0;i<p;i++){
    for (int j=0;j<K;j++){
      Cluster(i,j)=R::runif(Pmin(i,j),Pmax(i,j));
    }
  }

  return Cluster;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd COR(const Eigen::MatrixXd &X){

  const int p = X.cols();

  MatrixXd Cors(p,p);
  Cors = MatrixXd::Identity(p,p);

  for (int i = 1;i<p-1;i++){
    for (int j = i+1;j<p;j++){
      Cors(i,j)=X.col(i).transpose()*X.col(j);
      Cors(j,i)=Cors(i,j);
    }
  }

  return Cors;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
NumericMatrix CORYX(const NumericMatrix &Zs,
                    const NumericMatrix &Ys,
                    const NumericMatrix &Xs){

  Eigen::Map<Eigen::MatrixXd> Z = as<Eigen::Map<Eigen::MatrixXd> >(Zs);
  Eigen::Map<Eigen::MatrixXd> Y = as<Eigen::Map<Eigen::MatrixXd> >(Ys);
  Eigen::Map<Eigen::MatrixXd> X = as<Eigen::Map<Eigen::MatrixXd> >(Xs);

  const int q = Z.cols();
  const int p = Y.cols();
  const int r = X.cols();
  const int n = X.rows();
  const double dn = n;
  const int m = q+p+r;

  MatrixXd Cors  = MatrixXd::Zero(m,m);
  MatrixXd CorZY = MatrixXd::Zero(q,p);
  MatrixXd CorYX = MatrixXd::Zero(p,r);
  //Here I need to initiate Cors as a matrix of zeros

  for (int i = 0;i<q;i++){
    for (int j = 0;j<p;j++){
      CorZY(i,j)=Z.col(i).transpose()*Y.col(j);
    }
  }

  for (int i = 0;i<p;i++){
    for (int j = 0;j<r;j++){
      CorYX(i,j)=Y.col(i).transpose()*X.col(j);
    }
  }

  Cors.block(0,q,q,p) = CorZY/dn;
  Cors.block(q,0,p,q) = CorZY.transpose()/dn;
  Cors.block(q,q+p,p,r) = CorYX/dn;
  Cors.block(q+p,q,r,q) = CorYX.transpose()/dn;
  NumericMatrix Corsx(wrap(Cors));
  return Corsx;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
NumericMatrix CORYX2(const NumericMatrix &Zs,
                    const NumericMatrix &Ys,
                    const NumericMatrix &Xs){

  Eigen::Map<Eigen::MatrixXd> Z = as<Eigen::Map<Eigen::MatrixXd> >(Zs);
  Eigen::Map<Eigen::MatrixXd> Y = as<Eigen::Map<Eigen::MatrixXd> >(Ys);
  Eigen::Map<Eigen::MatrixXd> X = as<Eigen::Map<Eigen::MatrixXd> >(Xs);

  const int q = Z.cols();
  const int p = Y.cols();
  const int r = X.cols();
  const int n = X.rows();
  const double dn = n;
  const int m = q+p+r;

  MatrixXd Cors  = MatrixXd::Zero(m,m);
  MatrixXd CorZY = MatrixXd::Zero(q,p);
  MatrixXd CorYX = MatrixXd::Zero(p,r);
  //Here I need to initiate Cors as a matrix of zeros

  CorZY=Z.transpose()*Y;
  CorYX=Y.transpose()*X;

  Cors.block(0,q,q,p) = CorZY/dn;
  Cors.block(q,0,p,q) = CorZY.transpose()/dn;
  Cors.block(q,q+p,p,r) = CorYX/dn;
  Cors.block(q+p,q,r,q) = CorYX.transpose()/dn;

  NumericMatrix Corsx(wrap(Cors));
  return Corsx;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
NumericMatrix COR2(const NumericMatrix &Xs){

  Eigen::Map<Eigen::MatrixXd> X = as<Eigen::Map<Eigen::MatrixXd> >(Xs);
  const int p = X.cols();
  const double n = X.cols();

  MatrixXd Cors(p,p);
  Cors = MatrixXd::Identity(p,p);

  for (int i = 0;i<p-1;i++){
    for (int j = i+1;j<p;j++){
      //There is something wrong below
      Cors(i,j)=(X.col(i).transpose()*X.col(j)).sum()/n;
      Cors(j,i)=Cors(i,j);
    }
  }

  Rcpp::NumericMatrix Cor(wrap(Cors));
  return wrap(Cor);
}

// [[Rcpp::export]]
NumericMatrix matrixMAX(const NumericMatrix &A1,
                        const NumericMatrix &A2){

  const int p = A1.ncol();
  const int n = A1.nrow();

  NumericMatrix A(n,p);
  for (int i=0;i<n;i++){
    for(int j=0;j<p;j++){
      A(i,j)=max(NumericVector::create(A1(i,j),A2(i,j)));
    }
  }

  return A;
}

// [[Rcpp::export]]
NumericMatrix matrixMIN(const NumericMatrix &A1,
                        const NumericMatrix &A2){

  const int p = A1.ncol();
  const int n = A1.nrow();

  NumericMatrix A(n,p);
  for (int i=0;i<n;i++){
    for(int j=0;j<p;j++){
      A(i,j)=min(NumericVector::create(A1(i,j),A2(i,j)));
    }
  }

  return A;
}
