library(parallel)
library(microbenchmark)
library(profmem)

w = 5
LETTERS = c( 'a', 'c', 't', 'g' )
num_seq = 10
values <- 1:num_seq
nCpu <- 3
numWorkers <- nCpu

# Create sequence
genSeqs <- function(n = 10000, s=10) {
  do.call(paste0, replicate(s, sample(LETTERS, n, TRUE), FALSE))
}

SeqW <- genSeqs(num_seq, w)
SeqNextW <-  genSeqs(num_seq, w+1)

SeqNextWCut <- unlist(lapply(SeqNextW, function(x) substring(x, 2, 2*(w+1)+2-1)))

workerFunc<-function(n) { Find(length, list(grep(SeqW[n],SeqNextWCut,fixed=TRUE), 0))}

values <- 1:length(SeqW)

clusterProc <- function() {
  cl <- makeCluster(numWorkers)
  clusterExport(cl, "SeqW")
  clusterExport(cl, "SeqNextWCut")
  idx_2 <- parLapply(cl, 1:length(SeqW), function(n) {Find(length, list(grep(SeqW[n],SeqNextWCut,fixed=TRUE), 0))})
  stopCluster(cl)
}


Rcpp::cppFunction("std::vector<std::vector <int>> cpp_str_sort( std::vector< std::string > in_str, std::vector< std::string > out_str ) {
  #include <iostream>
  
  int len_in = in_str.size();
  int len_out = out_str.size();
  
  std::vector<std::vector<int>> indexes;

  for (int i = 0; i < len_in; i++){
    // std::vector <int> aux;
    for (int j = 0; j < len_out; j++){
      std::cout << j;
      //if (in_str[i].compare(out_str[j]) == 0){
      //  aux.push_back(j);
      //}
    }
    // indexes.push_back(aux);
  }
  return indexes;
}")

cpp_str_sort(SeqW, SeqNextWCut)


# microbenchmark(mclapply(values, workerFunc, mc.cores = numWorkers),
#                clusterProc(),
#                unit='s', times=10)