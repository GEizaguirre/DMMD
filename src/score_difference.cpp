#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector mean_diff_cpp(NumericVector bins1, NumericVector counts1, NumericVector bins2, NumericVector counts2) {
  
  size_t i;
  size_t num_bins1 = bins1.size();
  size_t num_bins2 = bins2.size();
  
  double sum_vals1 = 0;
  long int num_vals1 = 0;
  double mean_vals1;
  for ( i = 0; i < num_bins1; i++ ) {
    sum_vals1 += (bins1[i]*counts1[i]);
    num_vals1 += counts1[i];
  }
  mean_vals1 = sum_vals1/num_vals1;
  
  double sum_vals2 = 0;
  long int num_vals2 = 0;
  double mean_vals2;
  for ( i = 0; i < num_bins2; i++ ) {
    sum_vals2 += (bins2[i]*counts2[i]);
    num_vals2 += counts2[i];
  }
  mean_vals2 = sum_vals2/num_vals2;
  
  NumericVector result(2);
  result[0] = mean_vals1;
  result[1] = mean_vals2;
  
  return result;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
N = 1000000

# bins1 <- runif(N)
# counts1 <- sample.int(20, N, replace = TRUE)
# bins2 <- runif(N)
# counts2 <- sample.int(20, N, replace = TRUE)
# 
# res_cpp <- mean_diff_cpp(bins1, counts1, bins2, counts2)
# cdat_cpp <- data.frame(Set = c("Prone", "Resistant"), Score = c(res_cpp[0], res_cpp[1]))
# 
# mean_diff_R <- function(){
# 
#   p=unlist(sapply(1:length(bins1), function(i) rep(bins1[i],counts1[i])))
#   r=unlist(sapply(1:length(bins2), function(i) rep(bins2[i],counts2[i])))
# 
#   dat1 <- data.frame(Set = factor(rep(c("Prone"), each=length(p))), Score = p)
# 
#   dat2 <- data.frame(Set = factor(rep(c("Resistant"), each=length(r))), Score = r)
# 
#   dat <- rbind(dat1,dat2)
#   cdat <- ddply(dat, "Set", summarise, rating.mean=mean(Score))
# 
#   return(cdat)
# }
# 
# cdat_R <- mean_diff_R()

res <-mean_diff_cpp(pro_bin, pro_count, res_bin, res_count)

*/
