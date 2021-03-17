#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
List fuse_seqs_c( int motif_length, int grow_mode, std::vector<std::vector<int> > ind_fu,  std::vector<int> ind_fu_upd, 
                                            NumericVector fre_w, NumericVector sg_w, NumericMatrix fre_w_vec,
                                            NumericVector fre_w_next, NumericVector sg_w_next, NumericMatrix fre_w_vec_next) {
  
  int len_in = ind_fu_upd.size();
  int i1, fr, fr_next, len, idx, start_next;
  float sg, sg_next;
  
  for (int i = 0; i < len_in; i++){
    i1 = ind_fu_upd[i] - 1;
    fr = fre_w[i1];
    sg = sg_w[i1];
    
    len = ind_fu[i1].size();
    
    for (int k = 0; k < len; k++ ){
      idx = ind_fu[i1][k] - 1;
      fr_next = fre_w_next[idx];
      sg_next = sg_w_next[idx];

      sg_w_next[idx] = (fr*sg+fr_next*sg_next)/(fr+fr_next);
      switch(grow_mode){
      case 0: //C
        start_next = 1;
        break;
      case 1: //L
        start_next = 0;
        break;
      case 2: //R
        start_next = 2;
        break;
      }

      for (int j = start_next; j < (motif_length + start_next); j ++){
        fre_w_vec_next(idx, j) = fre_w_vec_next(idx, j) + fre_w_vec(i1, j-start_next);
        }

      fre_w_next[idx] = fr_next + fr;

      }
      
      
    }
  
  List ret;
  ret["sg_w_next"] = sg_w_next;
  ret["fre_w_next"] = fre_w_next;
  ret["fre_w_vec_next"] = fre_w_vec_next;
  return ret;
  
}


/*** R
w = 2
LETTERS = c( 'a', 'c', 't', 'g' )
num_seq = 10
values <- 1:num_seq
nCpu <- 3
numWorkers <- nCpu
'%ni%' = Negate('%in%') 

# Create sequence
genSeqs <- function(n = 10000, s=10) {
  do.call(paste0, replicate(s, sample(LETTERS, n, TRUE), FALSE))
}

SeqW <- genSeqs(num_seq, 2*w)
SeqNextW <-  genSeqs(num_seq, 2*w+2)
cut_point <- 2*w+1
SeqNextWCut <- unlist(lapply(SeqNextW, function(x) substring(x, 2, cut_point)))

print(SeqNextW)
print(SeqNextWCut)

# Get 
MetW = runif(num_seq)
print(MetW)
FreW = sample.int(20, num_seq, replace = TRUE)
print(FreW)
IndW = sample.int(20, num_seq, replace = TRUE)
print(IndW)
# FreWVec = FreVecW[[w]]

MetNextW = runif(num_seq)
FreNextW = sample.int(20, num_seq, replace = TRUE)
IndNextW = sample.int(20, num_seq, replace = TRUE)

IndFou_ <- list()
for (i in 1:num_seq){
  num_ind <- sample.int(10, 1, replace = TRUE)
  IndFou_[[i]] <- sort(sample.int(num_seq, num_ind, replace = FALSE))
}

IndFoUpd_ <- sort(sample.int(length(IndFou_), as.integer(2*length(IndFou_)/3), replace = FALSE))

FreWVec <- matrix(nrow=num_seq, ncol = 2*w)
for (i in 1:num_seq){
  FreWVec[i, ] <- sample.int(10, 2*w, replace = TRUE)
}

FreNextWVec <- matrix(nrow=num_seq, ncol = 2*w+2)
for (i in 1:num_seq){
  FreNextWVec[i, ] <- sample.int(10, 2*w+2, replace = TRUE)
}

grow_mode = 0

ret <- fuse_seqs_c(2*w, grow_mode, IndFou_, IndFoUpd_, FreW, MetW, FreWVec, FreNextW, MetNextW, FreNextWVec)

*/