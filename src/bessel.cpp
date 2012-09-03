

#include <ginac/ginac.h>

using namespace GiNaC;

#include <iostream>
#include <iomanip>

using namespace std;

#include "Rmath.h"


#include "bessel.h"




double eval_besselk(double x, double nu){
  return bessel_k(x,nu,1.);
}



#include <Rcpp.h>


RcppExport SEXP eval_besselk_R( SEXP x,SEXP nu){

  double xv = Rcpp::as<double>(x);
  double nuv = Rcpp::as<double>(nu);

  return Rcpp::wrap(eval_besselk(xv,nuv));
}


ex bessel_k_eval(const ex &x, const ex & nu){
  return bessel_k_fn(x,nu).hold();
}


ex bessel_k_evalf(const ex &x, const ex & nu){
  if( is_a<numeric>(x) && is_a<numeric>(nu) ) {
    return eval_besselk(ex_to<numeric>(x).to_double(),ex_to<numeric>(nu).to_double()) ;
  } else {
    return bessel_k_fn(x,nu).hold();
  }
}


ex bessel_k_eval_deri(const ex &x, const ex &nu,unsigned diffpar ){
  return - 0.5 * ( bessel_k_fn(x,nu+1) + bessel_k_fn(x,nu-1) );
}





REGISTER_FUNCTION(bessel_k_fn, 
 		  eval_func(bessel_k_eval).
 		  evalf_func(bessel_k_evalf).
  		  derivative_func(bessel_k_eval_deri)
		  );



