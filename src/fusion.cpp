#include <Rcpp.h> 
#include <omp.h>
using namespace Rcpp;


// [[Rcpp::export]]
void fuse_seqs_seq( int motif_length, int grow_mode, IntegerVector str_indexes, IntegerVector subindexes, IntegerVector ind_not_void, 
                    IntegerVector fre_w, NumericVector sg_w, IntegerMatrix fre_w_vec,
                    IntegerVector fre_w_next, NumericVector sg_w_next, IntegerMatrix fre_w_vec_next,
                    int num_cpu) {
  
  int len_next = ind_not_void.size();
  
  int i_start, i_end, k, i, ind_next, ind_current, j, start_next;
  int fr, fr_next;
  float sg, sg_next;
  
  for (i = 0; i < len_next; i++){
    
    ind_next = ind_not_void[i];
    i_start = subindexes[ind_next];
    i_end = subindexes[ind_next + 1];
    fr_next = fre_w_next[ind_next];
    sg_next = sg_w_next[ind_next];
    
    for (k = i_start; k < i_end; k++ ){
      
      ind_current = str_indexes[k];
      fr = fre_w[ind_current];
      sg = sg_w[ind_current];
      
      sg_w_next[ind_next] = (fr*sg+fr_next*sg_next)/(fr+fr_next);
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
      
      for (j = start_next; j < (motif_length + start_next); j ++){
        fre_w_vec_next(ind_next, j) = fre_w_vec_next(ind_next, j) + fre_w_vec(ind_current, j-start_next);
      }
      
      fre_w_next[ind_next] = fr_next + fr;
      
    }
    
  }
}

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
void fuse_seqs_openmp( int motif_length, int grow_mode, IntegerVector str_indexes, IntegerVector subindexes, IntegerVector ind_not_void, 
                       IntegerVector fre_w, NumericVector sg_w, IntegerMatrix fre_w_vec,
                       IntegerVector fre_w_next, NumericVector sg_w_next, IntegerMatrix fre_w_vec_next,
                       int num_cpu) {
  
  int len_next = ind_not_void.size();
  
  int i_start, i_end, k, i, ind_next, ind_current, j, start_next;
  int fr, fr_next;
  float sg, sg_next;
  
#pragma omp parallel num_threads(num_cpu)
#pragma omp for private(i_start, i_end, k, i, ind_next, ind_current, j, start_next, fr, fr_next, sg, sg_next)
  for (i = 0; i < len_next; i++){
    
    ind_next = ind_not_void[i];
    i_start = subindexes[ind_next];
    i_end = subindexes[ind_next + 1];
    fr_next = fre_w_next[ind_next];
    sg_next = sg_w_next[ind_next];
    
    for (k = i_start; k < i_end; k++ ){
      
      ind_current = str_indexes[k];
      fr = fre_w[ind_current];
      sg = sg_w[ind_current];
      
      sg_w_next[ind_next] = (fr*sg+fr_next*sg_next)/(fr+fr_next);
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
      
      for (j = start_next; j < (motif_length + start_next); j ++){
        fre_w_vec_next(ind_next, j) = fre_w_vec_next(ind_next, j) + fre_w_vec(ind_current, j-start_next);
      }
      
      fre_w_next[ind_next] = fr_next + fr;
      
    }
    
  }
  
}


// [[Rcpp::export]]
List find_strings_seq( StringVector in_str, StringVector out_str ) {
  
  int len_in = in_str.size();
  int len_out = out_str.size();
  
  int current_size = len_out*3;
  std::vector < int > str_indexes (current_size);
  std::vector < int > subindexes(len_out);
  std::vector < int > ind_not_void (len_out);
  
  int counter = 0;
  int prev_counter = counter;
  int not_void_counter = 0;

  for (size_t i = 0; i < len_out; i++){
    subindexes[i] = counter;
    
    for (size_t j = 0; j < len_in; j++){
      if (in_str(j) == out_str(i)){
        // resize index vector
        if (!( counter < current_size )){
          current_size = int(current_size + current_size / 2);
          str_indexes.resize(current_size);
          }
        str_indexes[counter] = j;
        counter += 1;
      }
    }
    if ( counter != prev_counter ){
      ind_not_void[not_void_counter++] = i;
      prev_counter = counter;
      }
    
  }
  
  str_indexes.resize(counter);
  ind_not_void.resize(not_void_counter);
  
  subindexes.push_back(str_indexes.size());
  
  std::vector < int > seqs_to_rm = str_indexes;
  std::sort(seqs_to_rm.begin(), seqs_to_rm.end());
  seqs_to_rm.erase(unique(seqs_to_rm.begin(),seqs_to_rm.end()),seqs_to_rm.end());
  
  for (int j = 0; j < seqs_to_rm.size(); j++ )
    seqs_to_rm[j]++;
  
  List ret;
  ret["str_indexes"] = str_indexes;
  ret["subindexes"] = subindexes;
  ret["ind_not_void"] = ind_not_void;
  ret["seqs_to_rm"] = seqs_to_rm;
  return ret;
}

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
List find_strings_par( StringVector in_str, StringVector out_str, int num_cpu ) {
  
  int len_in = in_str.size();
  int len_out = out_str.size();
  
  std::vector < int > str_indexes;
  std::vector < int > subindexes;
  std::vector < int > ind_not_void;
  std::vector < int > ids_per_thread (num_cpu);
  std::vector < int > ids_next_per_thread (num_cpu);
  
#pragma omp parallel num_threads(num_cpu)
{
  size_t i, j;
  int counter = 0;
  int next_seq_counter = 0;
  int prev_counter = counter;
  int current_size = int(len_out*3/num_cpu);
  int not_void_counter = 0;
  std::vector < int > str_indexes_priv (current_size);
  std::vector < int > ind_not_void_priv (int(len_out/num_cpu)+1);
  std::vector < int > subindexes_priv (int(len_out/num_cpu)+1);
  
#pragma omp for nowait schedule(static)
  for (i = 0; i < len_out; i++){
    
    subindexes_priv[next_seq_counter++] = counter;
    
    for (j = 0; j < len_in; j++){
      if (in_str(j) == out_str(i)){
        // resize index vector
        if (!( counter < current_size )){
          current_size = int(current_size + current_size / 2);
          str_indexes_priv.resize(current_size);
        }
        str_indexes_priv[counter++] = j;
      }
    }
    if ( counter != prev_counter ){
      ind_not_void_priv[not_void_counter++] = i;
      prev_counter = counter;
    }
  }
  
  str_indexes_priv.resize(counter);
  ind_not_void_priv.resize(not_void_counter);

  subindexes_priv.resize(next_seq_counter);

#pragma omp for schedule(static) ordered
  for(int i=0; i<omp_get_num_threads(); i++) {
#pragma omp ordered
    
    ids_next_per_thread[omp_get_thread_num()] = next_seq_counter;
    ids_per_thread[omp_get_thread_num()] = counter;
    ind_not_void.insert(ind_not_void.end(), ind_not_void_priv.begin(), ind_not_void_priv.end());
    str_indexes.insert(str_indexes.end(), str_indexes_priv.begin(), str_indexes_priv.end());
    subindexes.insert(subindexes.end(), subindexes_priv.begin(), subindexes_priv.end());
  }  
  
  }
  

  int accum_ids = ids_per_thread[0];
  int accum_ids_next = ids_next_per_thread[0];
  
  for (int idc = 1; idc < num_cpu; idc++ ){

    int end_id = accum_ids_next + ids_next_per_thread[idc];

    for (int ii = accum_ids_next; ii < end_id; ii++ ){
      subindexes[ii] += accum_ids;

    }

    accum_ids += ids_per_thread[idc];
    accum_ids_next = end_id;

  }
  
  subindexes.push_back(str_indexes.size());
  
  std::vector < int > seqs_to_rm = str_indexes;
  std::sort(seqs_to_rm.begin(), seqs_to_rm.end());
  seqs_to_rm.erase(unique(seqs_to_rm.begin(),seqs_to_rm.end()),seqs_to_rm.end());
  
  for (int j = 0; j < seqs_to_rm.size(); j++ )
    seqs_to_rm[j]++;
  
  List ret;
  ret["str_indexes"] = str_indexes;
  ret["subindexes"] = subindexes;
  ret["ind_not_void"] = ind_not_void;
  ret["seqs_to_rm"] = seqs_to_rm;
  return ret;
}

/*** R

w = 2
LETTERS = c( 'a', 'c', 't', 'g' )
num_seq = 10000
values <- 1:num_seq
nCpu <- 3
numWorkers <- nCpu
'%ni%' = Negate('%in%')

# Create sequence
genSeqs <- function(n = 10000, s=10) {
  do.call(paste0, replicate(s, sample(LETTERS, n, TRUE), FALSE))
}

parallel_exec <- function(){
  
  time1 <- Sys.time()
  SeqW <- genSeqs(num_seq, w)
  SeqNextW <-  genSeqs(num_seq, w+1)
  SeqNextWCut <- unlist(lapply(SeqNextW, function(x) substring(x, 2, 2*(w+1)+2-1)))
  
  SeqW <- genSeqs(num_seq, 2*w)
  SeqNextW <-  genSeqs(num_seq, 2*w+2)
  cut_point <- 2*w+1
  SeqNextWCut <- unlist(lapply(SeqNextW, function(x) substring(x, 2, cut_point)))
  
  # Get 
  MetW = runif(num_seq)
  FreW = sample.int(20, num_seq, replace = TRUE)
  IndW = sample.int(20, num_seq, replace = TRUE)
  
  MetNextW = runif(num_seq)
  FreNextW = sample.int(20, num_seq, replace = TRUE)
  IndNextW = sample.int(20, num_seq, replace = TRUE)
  
  FreWVec <- matrix(nrow=num_seq, ncol = 2*w)
  for (i in 1:num_seq){
    FreWVec[i, ] <- sample.int(10, 2*w, replace = TRUE)
  }
  
  FreNextWVec <- matrix(nrow=num_seq, ncol = 2*w+2)
  for (i in 1:num_seq){
    FreNextWVec[i, ] <- sample.int(10, 2*w+2, replace = TRUE)
  }
  time2 <- Sys.time()
  print(paste("data generated at", time2-time1, "s"))
  
  grow_mode = 0
  
  ind_data <- find_strings_seq(SeqW, SeqNextWCut)
  
  IndFu <- ind_data$str_indexes
  IndFu_sub <- ind_data$subindexes
  IndFu_not_void <- ind_data$ind_not_void
  seqs_to_rm <- ind_data$seqs_to_rm
  
  time3 <- Sys.time()
  print(paste("strings paired at", time3-time2, "s"))
  
  fuse_seqs_openmp(2*w, grow_mode,IndFu, IndFu_sub, IndFu_not_void,
                   FreW, MetW, FreWVec,
                   FreNextW, MetNextW, FreNextWVec,
                   nCpu)
  
  time4 <- Sys.time()
  print(paste("strings fusioned at", time4-time3, "s"))
}

sequential_exec <- function(){
  
  SeqW <- genSeqs(num_seq, 2*w)
  SeqNextW <-  genSeqs(num_seq, 2*w+2)
  cut_point <- 2*w+1
  SeqNextWCut <- unlist(lapply(SeqNextW, function(x) substring(x, 2, cut_point)))
  
  # Get 
  MetW = runif(num_seq)
  FreW = sample.int(20, num_seq, replace = TRUE)
  IndW = sample.int(20, num_seq, replace = TRUE)
  
  MetNextW = runif(num_seq)
  FreNextW = sample.int(20, num_seq, replace = TRUE)
  IndNextW = sample.int(20, num_seq, replace = TRUE)
  
  FreWVec <- matrix(nrow=num_seq, ncol = 2*w)
  for (i in 1:num_seq){
    FreWVec[i, ] <- sample.int(10, 2*w, replace = TRUE)
  }
  
  FreNextWVec <- matrix(nrow=num_seq, ncol = 2*w+2)
  for (i in 1:num_seq){
    FreNextWVec[i, ] <- sample.int(10, 2*w+2, replace = TRUE)
  }
  
  grow_mode = 0
  
  ind_data <- find_strings_seq(SeqW, SeqNextWCut)
  
  IndFu <- ind_data$str_indexes
  IndFu_sub <- ind_data$subindexes
  IndFu_not_void <- ind_data$ind_not_void
  seqs_to_rm <- ind_data$seqs_to_rm
  
}

SeqW <- genSeqs(num_seq, 2*w)
SeqNextW <-  genSeqs(num_seq, 2*w+2)
cut_point <- 2*w+1
SeqNextWCut <- unlist(lapply(SeqNextW, function(x) substring(x, 2, cut_point)))

data_par <- find_strings_par(SeqW, SeqNextWCut, nCpu)
data_seq <- find_strings_seq(SeqW, SeqNextWCut)

print(names(data_seq))

for (nm in names(data_seq)){
  print(paste("checking", nm))
  data1 <- data_par[[nm]]
  data2 <- data_seq[[nm]]
  if (length(data1) != length(data2)){
    print(paste("not equal size lists for ", nm))
    next
  }
  for (i in 1:length(data1)){
    if (data1[i] != data2[i]){
      print(paste("not equal lists ", nm))
      break
    }
  }
}

# library(microbenchmark)
# microbenchmark(
#   find_strings_par(SeqW, SeqNextWCut, nCpu),
#   find_strings_seq(SeqW, SeqNextWCut),
#   times = 10
# )


*/