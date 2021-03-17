#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<std::vector <int>> cpp_str_sort( StringVector in_str, StringVector out_str ) {
  
  int len_in = in_str.size();
  int len_out = out_str.size();
  
  std::vector<std::vector<int>> indexes;
  
  for (int i = 0; i < len_in; i++){
    std::vector <int> aux;
    for (int j = 0; j < len_out; j++){
      if (in_str(i) == out_str(j)){
        aux.push_back(j+1);
      }
    }
    if (aux.size()==0){
      aux.push_back(0);
    }
    indexes.push_back(aux);
  }
  return indexes;
}
/*std::vector<std::vector <int>> cpp_str_sort( std::vector< std::string > in_str, std::vector< std::string > out_str ) {
  #include <iostream>
  
  int len_in = in_str.size();
  int len_out = out_str.size();
  
  
  std::vector<std::vector<int>> indexes;
  
  for (int i = 0; i < len_in; i++){
    std::vector <int> aux;
    for (int j = 0; j < len_out; j++){
      if (in_str[i].compare(out_str[j]) == 0){
        aux.push_back(j+1);
      }
    }
    
    indexes.push_back(aux);
  }
  return indexes;
}*/


/*** R
w = 3
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

SeqW <- genSeqs(num_seq, w)
SeqNextW <-  genSeqs(num_seq, w+1)
SeqNextWCut <- unlist(lapply(SeqNextW, function(x) substring(x, 2, 2*(w+1)+2-1)))

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
# FreNextWVec=FreVecW[[w+1]]

IndFou <- cpp_str_sort(SeqW, SeqNextWCut)

IndFoUpd <- which(IndFou %ni% 0)

print(IndFou)
print(IndFoUpd)

# workerFunc<-function(n) { Find(length, list(grep(SeqW[n],SeqNextWCut,fixed=TRUE), 0))}
# lapply(1:num_seq, workerFunc)
*/