

#include <math.h>

#include <ginac/ginac.h>
#include <vector>

using namespace GiNaC;

#include "gtilde.partitions.h"
#include "bessel.h"


ex gtilde_1_new(const ex & x,int maxorder  = 20){
  std::vector<ex> exs;
  for( int i = 1 ; i <= maxorder ; i++){
    ex mulp = 4. * (double)getPartition(i) / ::sqrt( (double) i ) / x;
    exs.push_back( mulp * bessel_k_fn( ::sqrt((double)i) * x,1) );
  }

  return add(exs);
}



#include <Rcpp.h>


RcppExport SEXP eval_gtilde_n_R( SEXP x){

  double xv = Rcpp::as<double>(x);

  symbol sx;

  return Rcpp::wrap(ex_to<numeric>( gtilde_1_new(sx).subs(sx == xv).evalf()).to_double() ) ;
}



RcppExport SEXP eval_dgtilde_n_R( SEXP x){

  double xv = Rcpp::as<double>(x);

  symbol sx;

  return Rcpp::wrap(ex_to<numeric>( gtilde_1_new(sx).diff(sx).subs(sx == xv).evalf()).to_double() ) ;
}
