// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// fdr_c
float fdr_c(NumericMatrix ScoDatFraPro, NumericMatrix ScoDatFraRes, float lambda, int type_motif);
RcppExport SEXP _Motif_fdr_c(SEXP ScoDatFraProSEXP, SEXP ScoDatFraResSEXP, SEXP lambdaSEXP, SEXP type_motifSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ScoDatFraPro(ScoDatFraProSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ScoDatFraRes(ScoDatFraResSEXP);
    Rcpp::traits::input_parameter< float >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type type_motif(type_motifSEXP);
    rcpp_result_gen = Rcpp::wrap(fdr_c(ScoDatFraPro, ScoDatFraRes, lambda, type_motif));
    return rcpp_result_gen;
END_RCPP
}
// cpp_str_sort
std::vector<std::vector <int>> cpp_str_sort(StringVector in_str, StringVector out_str);
RcppExport SEXP _Motif_cpp_str_sort(SEXP in_strSEXP, SEXP out_strSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector >::type in_str(in_strSEXP);
    Rcpp::traits::input_parameter< StringVector >::type out_str(out_strSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_str_sort(in_str, out_str));
    return rcpp_result_gen;
END_RCPP
}
// fuse_seqs_c
List fuse_seqs_c(int motif_length, int grow_mode, std::vector<std::vector<int> > ind_fu, std::vector<int> ind_fu_upd, NumericVector fre_w, NumericVector sg_w, NumericMatrix fre_w_vec, NumericVector fre_w_next, NumericVector sg_w_next, NumericMatrix fre_w_vec_next);
RcppExport SEXP _Motif_fuse_seqs_c(SEXP motif_lengthSEXP, SEXP grow_modeSEXP, SEXP ind_fuSEXP, SEXP ind_fu_updSEXP, SEXP fre_wSEXP, SEXP sg_wSEXP, SEXP fre_w_vecSEXP, SEXP fre_w_nextSEXP, SEXP sg_w_nextSEXP, SEXP fre_w_vec_nextSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type motif_length(motif_lengthSEXP);
    Rcpp::traits::input_parameter< int >::type grow_mode(grow_modeSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<int> > >::type ind_fu(ind_fuSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type ind_fu_upd(ind_fu_updSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type fre_w(fre_wSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sg_w(sg_wSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type fre_w_vec(fre_w_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type fre_w_next(fre_w_nextSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sg_w_next(sg_w_nextSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type fre_w_vec_next(fre_w_vec_nextSEXP);
    rcpp_result_gen = Rcpp::wrap(fuse_seqs_c(motif_length, grow_mode, ind_fu, ind_fu_upd, fre_w, sg_w, fre_w_vec, fre_w_next, sg_w_next, fre_w_vec_next));
    return rcpp_result_gen;
END_RCPP
}
// scan_seqs_c
std::vector <float> scan_seqs_c(int nSeq, int LenMot, std::vector<std::vector<int>> NumSeq, NumericMatrix WeiLogPomElem, std::vector <float> WeiLogPWVElem);
RcppExport SEXP _Motif_scan_seqs_c(SEXP nSeqSEXP, SEXP LenMotSEXP, SEXP NumSeqSEXP, SEXP WeiLogPomElemSEXP, SEXP WeiLogPWVElemSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nSeq(nSeqSEXP);
    Rcpp::traits::input_parameter< int >::type LenMot(LenMotSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<int>> >::type NumSeq(NumSeqSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type WeiLogPomElem(WeiLogPomElemSEXP);
    Rcpp::traits::input_parameter< std::vector <float> >::type WeiLogPWVElem(WeiLogPWVElemSEXP);
    rcpp_result_gen = Rcpp::wrap(scan_seqs_c(nSeq, LenMot, NumSeq, WeiLogPomElem, WeiLogPWVElem));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP CooChr(SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP DissimilarityMatrix(SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP filtmdfile(SEXP, SEXP);
RcppExport SEXP readCooChrFile(SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP Reverse(SEXP, SEXP);
RcppExport SEXP scanPOMs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP SeqDic(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_Motif_fdr_c", (DL_FUNC) &_Motif_fdr_c, 4},
    {"_Motif_cpp_str_sort", (DL_FUNC) &_Motif_cpp_str_sort, 2},
    {"_Motif_fuse_seqs_c", (DL_FUNC) &_Motif_fuse_seqs_c, 10},
    {"_Motif_scan_seqs_c", (DL_FUNC) &_Motif_scan_seqs_c, 5},
    {"CooChr",              (DL_FUNC) &CooChr,              4},
    {"DissimilarityMatrix", (DL_FUNC) &DissimilarityMatrix, 4},
    {"filtmdfile",          (DL_FUNC) &filtmdfile,          2},
    {"readCooChrFile",      (DL_FUNC) &readCooChrFile,      4},
    {"Reverse",             (DL_FUNC) &Reverse,             2},
    {"scanPOMs",            (DL_FUNC) &scanPOMs,            8},
    {"SeqDic",              (DL_FUNC) &SeqDic,              8},
    {NULL, NULL, 0}
};

RcppExport void R_init_Motif(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
