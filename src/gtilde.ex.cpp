

#include <math.h>

#include <ginac/ginac.h>
#include <vector>

using namespace GiNaC;

#include "gtilde.partitions.h"
#include "bessel.h"

#include "gtilde.ex.h"


ex gtilde1_ex(const ex & x,int maxorder){
  std::vector<ex> exs;
  for( int i = 1 ; i <= maxorder ; i++){
    ex mulp = 4. * (double)getPartition(i) / ::sqrt( (double) i ) / x;
    exs.push_back( mulp * bessel_k_fn( ::sqrt((double)i) * x,1) );
  }

  return add(exs);
}



#include <Rcpp.h>


RcppExport SEXP eval_gtilde_ex_R( SEXP x){

  double xv = Rcpp::as<double>(x);

  symbol sx;

  return Rcpp::wrap(ex_to<numeric>( gtilde1_ex(sx).subs(sx == xv).evalf()).to_double() ) ;
}



RcppExport SEXP eval_dgtilde_ex_R( SEXP x){

  double xv = Rcpp::as<double>(x);

  symbol sx;

  return Rcpp::wrap(ex_to<numeric>( gtilde1_ex(sx).diff(sx).subs(sx == xv).evalf()).to_double() ) ;
}
