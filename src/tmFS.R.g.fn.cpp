
#include <math.h>
#include <ginac/ginac.h>

using namespace GiNaC;

#include <Rcpp.h>

using namespace Rcpp;

#include "tmFS.R.g.fn.h"

namespace chifit {


  /**
   * the \sigma(x) function : \sqrt{ 1 - \frac{4}{x} }
   */
  complex<double> R_sigma(complex<double> x){
    return sqrt(1.-4./x);
  }

  /**
   * the deri of \sigma(x) function
   */
  complex<double> R_dsigma(complex<double> x){
    return 2. / ( x * x * R_sigma( x ) );
  }

  /**
   * the 2nd deri of \sigma(x) function
   */
  complex<double> R_ddsigma(complex<double> x){
    complex<double> sigma = R_sigma(x);
    complex<double> dsigma = R_dsigma(x);
    return   - dsigma * dsigma*(  x * sigma + 1. / sigma ) ;
  }

  complex<double> R_g(complex<double> x){
    return R_sigma(x) * log( ( R_sigma(x) - 1. ) / (R_sigma(x)  + 1. ) ) + 2.;
  }


  complex<double> R_g_dlog(complex<double> x){
    - 1. / ( x * R_sigma(x) );
  }

  complex<double> R_dg(complex<double> x) {
    complex<double> sigma = R_sigma(x);
   return R_dsigma(x) * log( ( sigma - 1. ) / ( sigma + 1. ) ) - 1. / x;
  }

  complex<double> R_ddg(complex<double> x) {
    complex<double> sigma = R_sigma(x);
    return (
	    R_ddsigma(x) * log( ( sigma - 1. ) / ( sigma + 1. ) ) 
	    + 1. / ( x * x )*( 
			      1. - 2. / ( x * sigma * sigma )
			       )
	    );

  }


  RcppExport SEXP R_g_fn_R(SEXP x){
    return wrap(R_g(as<complex<double> >(x) ) );
  }

  RcppExport SEXP R_dg_fn_R(SEXP x){
    return wrap(R_dg(as<complex<double> >(x) ) );
  }

};
