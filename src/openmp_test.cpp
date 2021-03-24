#include <Rcpp.h>
#include <unistd.h>
#include <math.h>
#include <omp.h>
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export(c_bound_test_openmp)]]
NumericVector c_bound_test_openmp(NumericVector vin, int ncores)
{
  
int size_v = vin.size();
size_t i;
float aux;

NumericVector vout(size_v);
  
#pragma omp parallel num_threads(ncores)
#pragma omp for private(i, aux)
  for(i = 0; i < size_v; i++)
  { 
    aux = sqrt(pow(vin[i], 2)*2);
    vout[i] = aux + 3.14*aux;
  }
  
  return vout;
  
}

// [[Rcpp::export(c_bound_test_seq)]]
NumericVector c_bound_test_seq(NumericVector vin)
{
  
  int size_v = vin.size();
  float aux;
  NumericVector vout(size_v);

  for(size_t i = 0; i < size_v; i++)
  { 
    aux = sqrt(pow(vin[i], 2)*2);
    vout[i] = aux + 3.14*aux;
  }
  
  return vout;
  
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
ncores = 4
vin = runif(10000000)
library(microbenchmark)

microbenchmark(c_bound_test_seq(vin),
               c_bound_test_openmp(vin, 4),
               unit='s', times=10)
*/
