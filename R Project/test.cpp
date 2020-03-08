#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector gen_normal() {
  return rnorm(5000);
}