#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
float fdr_c(NumericMatrix ScoDatFraPro, NumericMatrix ScoDatFraRes, float lambda, int type_motif) {
  
  //TypeMotif = 0 (ForProne) = 1 (RevProne) = 2 (ForResis) = 3 (RevResis)
  int nrow_pro = ScoDatFraPro.nrow();
  int nrow_res = ScoDatFraRes.nrow();
  
  float sum = 0.0;
  int nelem = 0;
  float mean, std, aux, thr;

  // Calculate mean
  if ( type_motif < 2 ){
    for (int i = 0; i < nrow_pro; i++){
      sum += ScoDatFraPro(i, 0)*ScoDatFraPro(i, 1);
      nelem += ScoDatFraPro(i, 0);
    }
    mean = sum/nelem;
    sum = 0.0;
    for (int i = 0; i < nrow_pro; i++){
      aux = ScoDatFraPro(i, 1)-mean;
      sum += ScoDatFraPro(i, 0)*aux*aux;
    }
    std = sqrt(sum/(nelem-1));
  }
  else {
    for (int i = 0; i < nrow_res; i++){
      sum += ScoDatFraRes(i, 0)*ScoDatFraRes(i, 1);
      nelem += ScoDatFraRes(i, 0);
    }
    mean = sum/nelem;
    sum = 0.0;
    for (int i = 0; i < nrow_res; i++){
      aux = ScoDatFraRes(i, 1)-mean;
      sum += ScoDatFraRes(i, 0)*aux*aux;
    }
    std = sqrt(sum/(nelem-1));
  }
  
  thr = mean + lambda * std;
  
  int num_pro_above_thr = 0;
  int num_res_above_thr = 0;
  int intersection = (nrow_pro > nrow_res) ? nrow_res : nrow_pro;
  int smaller_set = (nrow_pro > nrow_res) ? 0 : 1;
  
  for (int i = 0; i < intersection; i++){
    if (ScoDatFraRes(i, 1) > thr){
      num_res_above_thr += ScoDatFraRes(i, 0);
    }
    if (ScoDatFraPro(i, 1) > thr){
      num_pro_above_thr += ScoDatFraPro(i, 0);
    }
  }
  
  if (smaller_set == 1){
    for (int i = intersection; i < nrow_res; i++){
      if (ScoDatFraRes(i, 1) > thr){
        num_res_above_thr += ScoDatFraRes(i, 0);
      }
    }
  }
  else {
    for (int i = intersection; i < nrow_pro; i++){
      if (ScoDatFraPro(i, 1) > thr){
        num_pro_above_thr += ScoDatFraPro(i, 0);
      }
    }
  }
  
  int tn = num_pro_above_thr + num_res_above_thr;
  float fdr_;
  
  if (tn == 0){
    fdr_ = 1.0;
  }
  else {
    if ( type_motif < 2 ) {
      fdr_ = float(num_res_above_thr) / float(tn);
    }
    else {
      fdr_ = float(num_pro_above_thr) / float(tn);
    }
    if (fdr_ == 0.0){
      fdr_ = 1;
    }
  }
  
  return fdr_;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
nelem1 <- 3000
nelem2 <- 3000
type_motif <- "RevResis"
type_motif_numeric <- switch(type_motif,
                           ForProne=0,
                           RevProne=1,
                           ForResis=2,
                           RevResis=3)

pro_reps <- sample.int(20, nelem1, replace = TRUE)
pro_vals <- runif(nelem1)
res_reps <- sample.int(20, nelem2, replace = TRUE)
res_vals <- runif(nelem2)

ScoDatFraPro <- cbind(pro_reps, pro_vals)
ScoDatFraRes <- cbind(res_reps, res_vals)
lambda <- 0.2


fdr_R <- function(ScoDatFraPro, ScoDatFraRes, lambda, TypeMotif){
  # For each <count,break> repeat the break count times.
  ScanProVec=unlist(sapply(1:length(ScoDatFraPro[,1]), function(i) rep(ScoDatFraPro[,2][i],ScoDatFraPro[,1][i])))
  ScanResVec=unlist(sapply(1:length(ScoDatFraRes[,1]), function(i) rep(ScoDatFraRes[,2][i],ScoDatFraRes[,1][i])))
  
  # Calculate the mean of the scores.
  Mea=switch(TypeMotif,
             ForProne=mean(ScanProVec),
             RevProne=mean(ScanProVec),
             ForResis=mean(ScanResVec),
             RevResis=mean(ScanResVec))
  #Log
  # line <- paste("Mean done for ", i, " of ", w, " from ", LenMotif)
  # write(line,file=Config$LogFile,append=TRUE)
  
  # Calculate the standard deviation of the scores.
  Std=switch(TypeMotif,
             ForProne=sd(ScanProVec),
             RevProne=sd(ScanProVec),
             ForResis=sd(ScanResVec),
             RevResis=sd(ScanResVec))
  
  #Log
  # line <- paste("sd done for ", i, " of ", w, " from ", LenMotif)
  # write(line,file=Config$LogFile,append=TRUE)
  
  # Calculate the threshold.
  Tr = Mea + lambda*Std
  
  # Select resistant and prone values.
  IndAboThrPro= which(ScanProVec>Tr)
  IndAboThrRes= which(ScanResVec>Tr)
  
  # Number of resistant and prone values.
  FP = length(IndAboThrPro)
  TP = length(IndAboThrRes)
  
  TN = FP +  TP
  
  # Calculate false discovery rate for this motif.
  fdr=switch(TypeMotif,
             ForProne = TP/TN,
             RevProne = TP/TN,
             ForResis = FP/TN,
             RevResis = FP/TN)
  
  if (is.na(fdr)==TRUE)  {
    fdr=1
  }
  else {
    if (fdr==0.0) fdr=1
  }
  
  return(fdr) 
}

# fdr_R_ <- fdr_R(ScoDatFraPro, ScoDatFraRes, lambda, type_motif)
# fdr_c_ <- fdr_c(ScoDatFraPro, ScoDatFraRes, lambda, type_motif_numeric)
# library(microbenchmark)
# microbenchmark(fdr_R(ScoDatFraPro, ScoDatFraRes, lambda, type_motif),
#                fdr_c(ScoDatFraPro, ScoDatFraRes, lambda, type_motif_numeric),
#                times = 10)

library(doMC)
# cl  <- makeCluster(4, 
#                    type = "SOCK")
# 
# src2 <- '
# float dumb_function(){
#   return 0.2;
# }
# '
# 
# clusterCall(cl, Rcpp::cppFunction, code=src2, env=environment())
# 
# 
# registerDoParallel(cl)
registerDoMC(4)

# clusterExport(cl, "fdr_c", envir = environment())
# clusterExport(cl, "dumb_function")

VecFdr=foreach(i=1:4,.combine=c) %dopar% {
  
  nelem1 <- 3000
  nelem2 <- 3000
  type_motif <- "RevResis"
  type_motif_numeric <- switch(type_motif,
                               ForProne=0,
                               RevProne=1,
                               ForResis=2,
                               RevResis=3)
  
  pro_reps <- sample.int(20, nelem1, replace = TRUE)
  pro_vals <- runif(nelem1)
  res_reps <- sample.int(20, nelem2, replace = TRUE)
  res_vals <- runif(nelem2)
  
  ScoDatFraPro <- cbind(pro_reps, pro_vals)
  ScoDatFraRes <- cbind(res_reps, res_vals)
  lambda <- 0.2
  
  type_motif_numeric
  
  # Get fdr value in cpp
  
  fdr <- fdr_c(ScoDatFraPro, ScoDatFraRes, lambda, type_motif_numeric)
  fdr

}

print(VecFdr)

*/
