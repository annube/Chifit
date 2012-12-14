

#include <ginac/ginac.h>
using namespace GiNaC;

#include <Rmath.h>
#include <math.h>

#include "utils.h"

#include "tmFS.B20.fn.h"

#include "gtilde.partitions.h"


#include <iostream>
#include <iomanip>




/*********************
 *
 * Evaluation of the B^2k functions
 *
 *
 *********************/


const int B2K_MAX_ORDER=40;

double eval_B2k_0(double x,int k){
  double sum=0.0;

  for( int i = 1 ; i <= B2K_MAX_ORDER ; i ++ ){
    double sx = sqrt( (double) i ) * x;
    double ds = (double) getPartition(i) *  bessel_k( sx , k + 1. , 1.) * pow( sx , (double) -(k+1) );
    sum += ds ; 
  }

  double stochastic_fact = 
    (double)( 2 * factorial(2*k))/ (double)( factorial(k) * exp2(k) );  // stochastik factor
  //  std::cout << "stoch. fact. " << stochastic_fact << std::endl;

  return 
    x 
    * stochastic_fact
    * sum;
}


/**
 * derivative of the above function w.r.t. x
 */

double eval_B2k_0_deri(double x,int k){
  double sum=0.0;

  for( int i = 1 ; i <= B2K_MAX_ORDER ; i ++ ){
    double sx = sqrt( (double) i ) * x;
    double ds = 
      - (double) getPartition(i)
      *  pow( sx , (double) -k ) 
      * (
	 ( 2. * (double) k + 1.)/sx * bessel_k(sx, k+1, 1.)
	 + bessel_k(sx, k , 1.)
	 );
    sum += ds ; 
  }

  double stochastic_fact = 
    (double)( 2 * factorial(2*k))/ (double)( factorial(k) * exp2(k) );  // stochastik factor
  //  std::cout << "stoch. fact. " << stochastic_fact << std::endl;

  return
    stochastic_fact
    * sum;
}




#include <Rcpp.h>
using namespace Rcpp;

RcppExport SEXP B2K_0_R(SEXP x,SEXP k){
  return wrap( eval_B2k_0 ( as<double>(x), as<int>(k) ) );
}

ex eval_B2k_0_fn(const ex & x,const ex &k){
  return B2k_0(x,k).hold();
}


ex evalf_B2k_0_fn(const ex & x,const ex &k){
  if( is_a<numeric>(x) &&  is_a<numeric>(k) ) {
    return eval_B2k_0(
		      ex_to<numeric>(x).to_double(),
		      MyExToInt(k)
		      ) ;
  } else {
    return B2k_0(x,k).hold();
  }
}

ex deri_B2k_0_fn(const ex & x,const ex &k,unsigned diffpar){
  if( diffpar == 0 )
    return B2k_0_deri(x,k).hold();
  else
    std::cerr << "implementation of derivative w.r.t. k not implemented.! " << std::endl;
  return NAN;
}


REGISTER_FUNCTION(B2k_0, 
 		  eval_func(eval_B2k_0_fn).
 		  evalf_func(evalf_B2k_0_fn).
		  derivative_func(deri_B2k_0_fn)
		  );


RcppExport SEXP dB2K_0_R(SEXP x,SEXP k){
  return wrap( eval_B2k_0_deri ( as<double>(x), as<int>(k) ) );
}


ex eval_B2k_0_deri_fn(const ex & x,const ex &k){
  return B2k_0_deri(x,k).hold();
}


ex evalf_B2k_0_deri_fn(const ex & x,const ex &k){
  if( is_a<numeric>(x) &&  is_a<numeric>(k) ) {
    return eval_B2k_0_deri(
		      ex_to<numeric>(x).to_double(),
		      MyExToInt(k)
		      ) ;
  } else {
    return B2k_0_deri(x,k).hold();
  }
}


REGISTER_FUNCTION(B2k_0_deri, 
 		  eval_func(eval_B2k_0_deri_fn).
 		  evalf_func(evalf_B2k_0_deri_fn)
		  );

