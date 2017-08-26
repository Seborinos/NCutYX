// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// NCutY3V1
double NCutY3V1(const NumericMatrix& Cys, const NumericMatrix& Cy2s, const NumericMatrix& Wys, const NumericMatrix& Wxs);
RcppExport SEXP NCutYX_NCutY3V1(SEXP CysSEXP, SEXP Cy2sSEXP, SEXP WysSEXP, SEXP WxsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Cys(CysSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Cy2s(Cy2sSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Wys(WysSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Wxs(WxsSEXP);
    rcpp_result_gen = Rcpp::wrap(NCutY3V1(Cys, Cy2s, Wys, Wxs));
    return rcpp_result_gen;
END_RCPP
}
// NCut
double NCut(const NumericMatrix& Cys, const NumericMatrix& Wys);
RcppExport SEXP NCutYX_NCut(SEXP CysSEXP, SEXP WysSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Cys(CysSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Wys(WysSEXP);
    rcpp_result_gen = Rcpp::wrap(NCut(Cys, Wys));
    return rcpp_result_gen;
END_RCPP
}
// WNCut
double WNCut(const NumericMatrix& Cys, const NumericMatrix& Cy2s, const NumericMatrix& Wys);
RcppExport SEXP NCutYX_WNCut(SEXP CysSEXP, SEXP Cy2sSEXP, SEXP WysSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Cys(CysSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Cy2s(Cy2sSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Wys(WysSEXP);
    rcpp_result_gen = Rcpp::wrap(WNCut(Cys, Cy2s, Wys));
    return rcpp_result_gen;
END_RCPP
}
// WNCut2
double WNCut2(const NumericMatrix& Cys, const NumericMatrix& Cy2s, const NumericMatrix& Dys, const NumericMatrix& Wys);
RcppExport SEXP NCutYX_WNCut2(SEXP CysSEXP, SEXP Cy2sSEXP, SEXP DysSEXP, SEXP WysSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Cys(CysSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Cy2s(Cy2sSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Dys(DysSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Wys(WysSEXP);
    rcpp_result_gen = Rcpp::wrap(WNCut2(Cys, Cy2s, Dys, Wys));
    return rcpp_result_gen;
END_RCPP
}
// WNCut3
double WNCut3(const NumericMatrix& Cys, const NumericMatrix& Cy2s, const NumericMatrix& Wys);
RcppExport SEXP NCutYX_WNCut3(SEXP CysSEXP, SEXP Cy2sSEXP, SEXP WysSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Cys(CysSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Cy2s(Cy2sSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Wys(WysSEXP);
    rcpp_result_gen = Rcpp::wrap(WNCut3(Cys, Cy2s, Wys));
    return rcpp_result_gen;
END_RCPP
}
// WNCut4
double WNCut4(const NumericMatrix& Cys, const NumericMatrix& Cy2s, const NumericMatrix& Dys, const NumericMatrix& Wys);
RcppExport SEXP NCutYX_WNCut4(SEXP CysSEXP, SEXP Cy2sSEXP, SEXP DysSEXP, SEXP WysSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Cys(CysSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Cy2s(Cy2sSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Dys(DysSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Wys(WysSEXP);
    rcpp_result_gen = Rcpp::wrap(WNCut4(Cys, Cy2s, Dys, Wys));
    return rcpp_result_gen;
END_RCPP
}
// WNCut5
double WNCut5(const NumericMatrix& Cys, const NumericMatrix& Dys, const NumericMatrix& Wys);
RcppExport SEXP NCutYX_WNCut5(SEXP CysSEXP, SEXP DysSEXP, SEXP WysSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Cys(CysSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Dys(DysSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Wys(WysSEXP);
    rcpp_result_gen = Rcpp::wrap(WNCut5(Cys, Dys, Wys));
    return rcpp_result_gen;
END_RCPP
}
// WNCut6
double WNCut6(const NumericMatrix& Cys, const NumericMatrix& Cy2s, const NumericMatrix& Wys);
RcppExport SEXP NCutYX_WNCut6(SEXP CysSEXP, SEXP Cy2sSEXP, SEXP WysSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Cys(CysSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Cy2s(Cy2sSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Wys(WysSEXP);
    rcpp_result_gen = Rcpp::wrap(WNCut6(Cys, Cy2s, Wys));
    return rcpp_result_gen;
END_RCPP
}
// WNCut7
double WNCut7(const NumericMatrix& Cys, const NumericMatrix& Cy2s, const NumericMatrix& Dys);
RcppExport SEXP NCutYX_WNCut7(SEXP CysSEXP, SEXP Cy2sSEXP, SEXP DysSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Cys(CysSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Cy2s(Cy2sSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Dys(DysSEXP);
    rcpp_result_gen = Rcpp::wrap(WNCut7(Cys, Cy2s, Dys));
    return rcpp_result_gen;
END_RCPP
}
// WNCut8
double WNCut8(const NumericMatrix& Cys, const NumericMatrix& Cy2s, const NumericMatrix& Wys);
RcppExport SEXP NCutYX_WNCut8(SEXP CysSEXP, SEXP Cy2sSEXP, SEXP WysSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Cys(CysSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Cy2s(Cy2sSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Wys(WysSEXP);
    rcpp_result_gen = Rcpp::wrap(WNCut8(Cys, Cy2s, Wys));
    return rcpp_result_gen;
END_RCPP
}
// NCutLayer3V1
double NCutLayer3V1(const NumericMatrix& Cys, const NumericMatrix& Cy2s, const NumericMatrix& Wzs, const NumericMatrix& Wys, const NumericMatrix& Wxs, const NumericMatrix& Wzyxs);
RcppExport SEXP NCutYX_NCutLayer3V1(SEXP CysSEXP, SEXP Cy2sSEXP, SEXP WzsSEXP, SEXP WysSEXP, SEXP WxsSEXP, SEXP WzyxsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Cys(CysSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Cy2s(Cy2sSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Wzs(WzsSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Wys(WysSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Wxs(WxsSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Wzyxs(WzyxsSEXP);
    rcpp_result_gen = Rcpp::wrap(NCutLayer3V1(Cys, Cy2s, Wzs, Wys, Wxs, Wzyxs));
    return rcpp_result_gen;
END_RCPP
}
// Penal
NumericMatrix Penal(const NumericMatrix& Cys);
RcppExport SEXP NCutYX_Penal(SEXP CysSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Cys(CysSEXP);
    rcpp_result_gen = Rcpp::wrap(Penal(Cys));
    return rcpp_result_gen;
END_RCPP
}
// Ranking
double Ranking(const NumericMatrix& C);
RcppExport SEXP NCutYX_Ranking(SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(Ranking(C));
    return rcpp_result_gen;
END_RCPP
}
// Ranking2
double Ranking2(const NumericMatrix& C);
RcppExport SEXP NCutYX_Ranking2(SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(Ranking2(C));
    return rcpp_result_gen;
END_RCPP
}
// Ranking3
double Ranking3(const NumericMatrix& C);
RcppExport SEXP NCutYX_Ranking3(SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(Ranking3(C));
    return rcpp_result_gen;
END_RCPP
}
// Ranking4
double Ranking4(const NumericMatrix& C);
RcppExport SEXP NCutYX_Ranking4(SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(Ranking4(C));
    return rcpp_result_gen;
END_RCPP
}
// Ranking5
double Ranking5(const NumericMatrix& C);
RcppExport SEXP NCutYX_Ranking5(SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(Ranking5(C));
    return rcpp_result_gen;
END_RCPP
}
// Ranking6
double Ranking6(const NumericMatrix& C, const double alpha);
RcppExport SEXP NCutYX_Ranking6(SEXP CSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(Ranking6(C, alpha));
    return rcpp_result_gen;
END_RCPP
}
// oneMultinomCalt
IntegerVector oneMultinomCalt(NumericVector probs);
RcppExport SEXP NCutYX_oneMultinomCalt(SEXP probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type probs(probsSEXP);
    rcpp_result_gen = Rcpp::wrap(oneMultinomCalt(probs));
    return rcpp_result_gen;
END_RCPP
}
// RandomMatrix
IntegerMatrix RandomMatrix(const int& p, const int& K, const NumericMatrix& P);
RcppExport SEXP NCutYX_RandomMatrix(SEXP pSEXP, SEXP KSEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(RandomMatrix(p, K, P));
    return rcpp_result_gen;
END_RCPP
}
// RandomUnifMatrix
NumericMatrix RandomUnifMatrix(const int& p, const int& K, const NumericMatrix& Pmin, const NumericMatrix& Pmax);
RcppExport SEXP NCutYX_RandomUnifMatrix(SEXP pSEXP, SEXP KSEXP, SEXP PminSEXP, SEXP PmaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Pmin(PminSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Pmax(PmaxSEXP);
    rcpp_result_gen = Rcpp::wrap(RandomUnifMatrix(p, K, Pmin, Pmax));
    return rcpp_result_gen;
END_RCPP
}
// COR
Eigen::MatrixXd COR(const Eigen::MatrixXd& X);
RcppExport SEXP NCutYX_COR(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(COR(X));
    return rcpp_result_gen;
END_RCPP
}
// COR2
NumericMatrix COR2(const NumericMatrix& Xs);
RcppExport SEXP NCutYX_COR2(SEXP XsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Xs(XsSEXP);
    rcpp_result_gen = Rcpp::wrap(COR2(Xs));
    return rcpp_result_gen;
END_RCPP
}
// matrixMAX
NumericMatrix matrixMAX(const NumericMatrix& A1, const NumericMatrix& A2);
RcppExport SEXP NCutYX_matrixMAX(SEXP A1SEXP, SEXP A2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type A1(A1SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type A2(A2SEXP);
    rcpp_result_gen = Rcpp::wrap(matrixMAX(A1, A2));
    return rcpp_result_gen;
END_RCPP
}
// matrixMIN
NumericMatrix matrixMIN(const NumericMatrix& A1, const NumericMatrix& A2);
RcppExport SEXP NCutYX_matrixMIN(SEXP A1SEXP, SEXP A2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type A1(A1SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type A2(A2SEXP);
    rcpp_result_gen = Rcpp::wrap(matrixMIN(A1, A2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"NCutYX_NCutY3V1", (DL_FUNC) &NCutYX_NCutY3V1, 4},
    {"NCutYX_NCut", (DL_FUNC) &NCutYX_NCut, 2},
    {"NCutYX_WNCut", (DL_FUNC) &NCutYX_WNCut, 3},
    {"NCutYX_WNCut2", (DL_FUNC) &NCutYX_WNCut2, 4},
    {"NCutYX_WNCut3", (DL_FUNC) &NCutYX_WNCut3, 3},
    {"NCutYX_WNCut4", (DL_FUNC) &NCutYX_WNCut4, 4},
    {"NCutYX_WNCut5", (DL_FUNC) &NCutYX_WNCut5, 3},
    {"NCutYX_WNCut6", (DL_FUNC) &NCutYX_WNCut6, 3},
    {"NCutYX_WNCut7", (DL_FUNC) &NCutYX_WNCut7, 3},
    {"NCutYX_WNCut8", (DL_FUNC) &NCutYX_WNCut8, 3},
    {"NCutYX_NCutLayer3V1", (DL_FUNC) &NCutYX_NCutLayer3V1, 6},
    {"NCutYX_Penal", (DL_FUNC) &NCutYX_Penal, 1},
    {"NCutYX_Ranking", (DL_FUNC) &NCutYX_Ranking, 1},
    {"NCutYX_Ranking2", (DL_FUNC) &NCutYX_Ranking2, 1},
    {"NCutYX_Ranking3", (DL_FUNC) &NCutYX_Ranking3, 1},
    {"NCutYX_Ranking4", (DL_FUNC) &NCutYX_Ranking4, 1},
    {"NCutYX_Ranking5", (DL_FUNC) &NCutYX_Ranking5, 1},
    {"NCutYX_Ranking6", (DL_FUNC) &NCutYX_Ranking6, 2},
    {"NCutYX_oneMultinomCalt", (DL_FUNC) &NCutYX_oneMultinomCalt, 1},
    {"NCutYX_RandomMatrix", (DL_FUNC) &NCutYX_RandomMatrix, 3},
    {"NCutYX_RandomUnifMatrix", (DL_FUNC) &NCutYX_RandomUnifMatrix, 4},
    {"NCutYX_COR", (DL_FUNC) &NCutYX_COR, 1},
    {"NCutYX_COR2", (DL_FUNC) &NCutYX_COR2, 1},
    {"NCutYX_matrixMAX", (DL_FUNC) &NCutYX_matrixMAX, 2},
    {"NCutYX_matrixMIN", (DL_FUNC) &NCutYX_matrixMIN, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_NCutYX(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
