#include <RcppEigen.h>
#include <algorithm>  // std::sort std::min

using namespace Rcpp;
using namespace Eigen;

//This function calculates the NCutYX objective function.

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double NCutY3V1(const NumericMatrix &Cys,
                const NumericMatrix &Cy2s,
                const NumericMatrix &Wys,
                const NumericMatrix &Wxs){

  Eigen::Map<Eigen::MatrixXd> Wy  = as<Eigen::Map<Eigen::MatrixXd> >(Wys);
  Eigen::Map<Eigen::MatrixXd> Wx  = as<Eigen::Map<Eigen::MatrixXd> >(Wxs);
  Eigen::Map<Eigen::MatrixXd> Cy  = as<Eigen::Map<Eigen::MatrixXd> >(Cys);
  Eigen::Map<Eigen::MatrixXd> Cy2 = as<Eigen::Map<Eigen::MatrixXd> >(Cy2s);
  const int K = Cy.cols();//number of groups
  //These vectors wil contain the number of vertices in from X and Y
  //For the K=2 case, cut1=cut2
  VectorXd Cuty(K);//Numerator of NCut
  VectorXd Volx(K);//Denoinator of NCut
  double J;
  for(int i = 0;i<K;i++){
    Cuty(i) = Cy.col(i).transpose()*Wy*Cy2.col(i);//This is new
    Volx(i) = Cy.col(i).transpose()*Wx*Cy.col(i);
  }

  VectorXd Res = Cuty.array()/Volx.array();
  J = Res.sum();

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
    Cuty(i) = (Cy.col(i).transpose()*Wy*Cy2.col(i)).sum();
    Volx(i) = (Cy.col(i).transpose()*Wy*Cy.col(i)).sum();
  }

  VectorXd Res=Cuty.array()/Volx.array();
  J=Res.sum();
  return J;
}

//This function calculates a NCut for a three layered graph.
//It has as input the cluster membership, and the respective
//weigth matrices.

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double NCutLayer3V1(const NumericMatrix &Cys,
                    const NumericMatrix &Cy2s,
                    const NumericMatrix &Wzs,
                    const NumericMatrix &Wys,
                    const NumericMatrix &Wxs,
                    const NumericMatrix &Wzyxs){

  Eigen::Map<Eigen::MatrixXd> Wz   = as<Eigen::Map<Eigen::MatrixXd> >(Wzs);
  Eigen::Map<Eigen::MatrixXd> Wy   = as<Eigen::Map<Eigen::MatrixXd> >(Wys);
  Eigen::Map<Eigen::MatrixXd> Wx   = as<Eigen::Map<Eigen::MatrixXd> >(Wxs);
  Eigen::Map<Eigen::MatrixXd> Wzyx = as<Eigen::Map<Eigen::MatrixXd> >(Wzyxs);
  Eigen::Map<Eigen::MatrixXd> Cy   = as<Eigen::Map<Eigen::MatrixXd> >(Cys);
  Eigen::Map<Eigen::MatrixXd> Cy2  = as<Eigen::Map<Eigen::MatrixXd> >(Cy2s);

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
  return J;
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

  for (int i=0;i<p;i++){
        Cluster(i,_)=oneMultinomCalt(P(i,_));
  }

  return Cluster;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
NumericMatrix RandomMatrix2(const int &p,
                            const int &K,
                            const NumericMatrix &P){
  NumericMatrix Cluster(p,K);

  for (int i=0;i<p;i++){
    Cluster(i,_)=oneMultinomCalt(P(i,_));
  }

  return Cluster;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List samplingncut(const NumericMatrix &W,
                  const NumericMatrix &Prob,
                  const int &p,
                  const int &K,
                  const int &N){
  NumericVector l(N);
  List C(N);
  for (int k = 0; k < N; k++){
    // draw from Prob a Cluster and calculate its loss
    C[k] = RandomMatrix2(p,K,Prob);
    l[k] = NCut(C[k],W);
  }
  return List::create(Named("Clusters") = C,
                      Named("loss")     = l);
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double cutoff(NumericVector &loss,
              const int &q0){
  std::sort(loss.begin(), loss.end());
  double qloss = loss(q0);
  return qloss;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
NumericMatrix ProbAve(List &Cs,
                      IntegerVector &Ind,
                      const int &p,
                      const int &K){
  NumericMatrix Prob(p,K);
  Rcpp::Function Reduce("Reduce");
  Prob = Reduce('+',Cs[Ind]);
  double N2 = Ind.length();
  double w = 1.0/N2;
  Prob = Prob*w;
  return Prob;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
IntegerVector Indexing(NumericVector &loss,
                       const double &qloss){
  // create an index Ind
  Rcpp::Function which("which");
  IntegerVector Ind = which(loss<=qloss);
  return Ind;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List ncutcem(const NumericMatrix &W,
                      const int &p,
                      const int &K,
                      const int &N,
                      const int &B,
                      const int &q0,
                      const double &p0){
  NumericMatrix Prob(p,K);
  NumericVector loss(N);
  NumericVector Quant(B);
  NumericVector loss2(N);
  List          Cs(N);
  // start some intial distribution for the variables
  std::fill(Prob.begin(), Prob.end(), p0);
  for (int i = 0; i < B; i++){
    List sample = samplingncut(W, Prob, p, K, N);
    Cs          = sample["Clusters"];
    loss        = sample["loss"];
    loss2       = clone(loss);
    // calculate the quantile and create an index
    double qloss      = cutoff(loss2, q0);
    Quant(i)          = qloss;
    IntegerVector Ind = Indexing(loss, qloss);
    Ind               = Ind - 1;
    Prob              = ProbAve(Cs, Ind, p, K);
  }
  return List::create(Named("quantile") = Quant,
                      Named("clusters") = Prob,
                      Named("ncut")     = NCut(Prob, W));
}

