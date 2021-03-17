#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector <float> scan_seqs_c(int nSeq, int LenMot, std::vector<std::vector<int>> NumSeq, NumericMatrix WeiLogPomElem, std::vector <float> WeiLogPWVElem) {
  
  std::vector <float> ScoRes(nSeq);
  float sco;
  int idx_aux;
  
  for (int i = 0; i < nSeq; i++){
    sco = 0;
    for (int j = 0; j < LenMot; j++ ){
      // Get the nucleotide in the current position of the sequence.
      // Get the score for the nucleotide in the motif.
      // Get the normalization value for the position.
      sco = sco + WeiLogPomElem(NumSeq[i][j], j) - WeiLogPWVElem[j];
    }
    ScoRes[i] = sco;
  }
  
  return ScoRes;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
NumSeq <- list(c(0,3,2,0,1), c(2,3,1,3,0), c(0,0,1,1,0))
WeiLogPWVElem <- runif(5)
WeiLogPomElem<-matrix(runif(20), nrow=4, ncol=5)
print(WeiLogPomElem)
scan_seqs(3, 5, NumSeq, WeiLogPomElem, WeiLogPWVElem)
*/
