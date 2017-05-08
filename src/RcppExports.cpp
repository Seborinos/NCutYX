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
// WNCut
double WNCut(const NumericMatrix& Cys, const NumericMatrix& Cy2s, const NumericMatrix& DWys, const NumericMatrix& Wys);
RcppExport SEXP NCutYX_WNCut(SEXP CysSEXP, SEXP Cy2sSEXP, SEXP DWysSEXP, SEXP WysSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Cys(CysSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Cy2s(Cy2sSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type DWys(DWysSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Wys(WysSEXP);
    rcpp_result_gen = Rcpp::wrap(WNCut(Cys, Cy2s, DWys, Wys));
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
