

/**
 * this source file will contain the "R_0^k" function
 * appearing in the colangelo-wenger-wu paper
 */



#include <complex>

#include <iostream>
#include <iomanip>

using namespace std;


#include <ginac/ginac.h>

using namespace GiNaC;



#include <gsl/gsl_integration.h>

#include <Rmath.h>


#include "utils.h"
#include "tmFS.R.g.fn.h"
#include "tmFS.R.fn.h"
#include "gtilde.partitions.h"






/*********************
 *
 * Evaluation of the R_0^k functions
 *
 *
 *********************/

typedef struct R_Integrand_Params_ {
  int k,n;
  double x,r;


} R_Integrand_Params;

double R_integrand(double y, void *params){
  //  static int fnEvalCount=0;
  R_Integrand_Params *sp=(R_Integrand_Params *)params;
  //  ++fnEvalCount;
  //  cout << fnEvalCount << endl;


  if( sp->k % 2 == 0 ){
    return 
      pow(y,sp->k)
      * exp( -sp->x * sqrt( sp->n* ( 1. + y * y ))  )
      *
      chifit::R_g( complex<double>(2.,2.*sp->r*y) ).real();
  } else {
    return 
      pow(y,sp->k)
      * exp( -sp->x * sqrt( sp->n* (1. + y * y) )  )
      *
      chifit::R_g( complex<double>(2.,2.*sp->r*y) ).imag();
  }
    
}

double R_integrand_dx(double y, void *params){
  //  static int fnEvalCount=0;
  R_Integrand_Params *sp=(R_Integrand_Params *)params;
  //  ++fnEvalCount;
  //  cout << fnEvalCount << endl;


  if( sp->k % 2 == 0 ){
    return 
      - pow(y,sp->k) * sqrt( sp->n*(1. + y * y) )
      * exp( -sp->x * sqrt( sp->n*(1. + y * y ))  )
      *
      chifit::R_g( complex<double>(2.,2.*sp->r*y) ).real();
  } else {
    return 
      - pow(y,sp->k) * sqrt( sp->n* ( 1. + y * y ) )
      * exp( -sp->x * sqrt( sp->n* ( 1. + y * y )  ) )
      *
      chifit::R_g( complex<double>(2.,2.*sp->r*y) ).imag();
  }
    
}

double R_integrand_dr(double y, void *params){
  //  static int fnEvalCount=0;
  R_Integrand_Params *sp=(R_Integrand_Params *)params;
  //  ++fnEvalCount;
  //  cout << fnEvalCount << endl;


  if( sp->k % 2 == 0 ){
    return 
      pow(y,sp->k) 
      * exp( -sp->x * sqrt( sp->n* ( 1. + y * y )  ) ) 
      *
      ( chifit::R_dg( complex<double>(2.,2.*sp->r*y) ) * complex<double>(0,2*y) ).real();
  } else {
    return 
      pow(y,sp->k)
      * exp( -sp->x * sqrt( sp->n* ( 1. + y * y ) ) )
      *
      ( chifit::R_dg( complex<double>(2.,2.*sp->r*y) ) * complex<double>(0,2*y) ).imag();
  }
    
}


//this is the implementation of the R function devided by r_i^(k+1) with the r_i substituted out of the exponent
// to get the correct R function you have to multiply this function with r_i^(k+1) !!!


double eval_tmFS_R(double x,int k,double r,IntegrandFN integrand=&R_integrand){

//   cout << "k = " << (k) << endl
//        << "x = " << (x) << endl
//        << "r = " << (r) << endl ;

  int WS_SIZE=1000;
  const int MAX_ORDER=40;

  R_Integrand_Params rip;
  rip.k=(k);
  rip.x=(x);
  rip.r=(r);

  static  gsl_integration_workspace *w=NULL;
  if( w == NULL)
    w = gsl_integration_workspace_alloc(WS_SIZE);

  gsl_function F;
  F.function = integrand;
  F.params = &rip;


  double intResult,intAbsErr;

  double sum=0;

  for( int n = MAX_ORDER ; n >= 1 ; n--){
    rip.x = x ;
    rip.n = n ;
    gsl_integration_qagi(&F,0,1.e-6,WS_SIZE,w,&intResult,&intAbsErr);
//     cout << "Integration result : " << intResult << endl;
//     cout << "Abs Error          : " << intAbsErr << endl;
//     cout << " intervals         : " << w->size << endl;

//    double reldiff = intResult/sum;
//    cout << "ds/s = " << intResult/sum << endl;

    sum += (double) getPartition(n)/sqrt((double)n) * intResult ;

  }


//  gsl_integration_workspace_free(w);

  return (sum);

}



#include <Rcpp.h>

using namespace Rcpp;

RcppExport SEXP tmFS_R_fn_R(SEXP x,SEXP k, SEXP r){
  return wrap( eval_tmFS_R(as<double>(x),as<int>(k), as<double>(r) ) );
}

RcppExport SEXP tmFS_R_dx_fn_R(SEXP x,SEXP k, SEXP r){
  return wrap( eval_tmFS_R(as<double>(x),as<int>(k), as<double>(r),&R_integrand_dx ) );
}



ex eval_tmFS_R_fn(const ex &x,const ex & k , const ex & r){
  if( is_a<numeric>(k) ) 
    return tmFS_R(x,ex_to<numeric>(k).to_int(),r).hold();
  else 
    return tmFS_R(x,k,r).hold();
}

ex eval_tmFS_R_fn_deri(const ex &x,const ex & k , const ex & r,unsigned diffpar){
  if(diffpar == 0 )
    return tmFS_R_dx(x,k,r);
  else if( diffpar == 2 ){
    cout << "Warning derivative of R w.r.t. r not implemented so far " << endl;
  } else {
    cerr << "Error : R can only be derived w.r.t. x and r " << endl;
  }
  return NAN;
}


ex evalf_tmFS_R_fn(const ex &x,const ex & k , const ex & r){
  if( is_a<numeric>(x) &&  is_a<numeric>(k) && is_a<numeric>(r) ) {

    return eval_tmFS_R(
		      ex_to<numeric>(x).to_double(),
		      MyExToInt( k ),
		      ex_to<numeric>(r).to_double()
		      ) ;
  } else {
    return tmFS_R(x,k,r).hold();
  }
}


RcppExport SEXP tmFS_R_fn_R_ginac(SEXP x,SEXP k, SEXP r){
  ex fnv=tmFS_R(as<double>(x),as<int>(k),as<double>(r));
  return wrap( ex_to<numeric>(fnv.evalf()).to_double() );
}



ex eval_tmFS_R_dx_fn(const ex &x,const ex & k , const ex & r){
  if( is_a<numeric>(k) ) 
    return tmFS_R(x,ex_to<numeric>(k).to_int(),r).hold();
  else 
    return tmFS_R(x,k,r).hold();
}


ex evalf_tmFS_R_dx_fn(const ex &x,const ex & k , const ex & r){
  if( is_a<numeric>(x) &&  is_a<numeric>(k) && is_a<numeric>(r) ) {

    return eval_tmFS_R(
		      ex_to<numeric>(x).to_double(),
		      MyExToInt( k ),
		      ex_to<numeric>(r).to_double(),
		      &R_integrand_dx
		      ) ;
  } else {
    return tmFS_R(x,k,r).hold();
  }
}




REGISTER_FUNCTION(tmFS_R,
		  eval_func(eval_tmFS_R_fn).
		  evalf_func(evalf_tmFS_R_fn).
		  derivative_func(eval_tmFS_R_fn_deri)
		  );

REGISTER_FUNCTION(tmFS_R_dx,
		  eval_func(eval_tmFS_R_dx_fn).
		  evalf_func(evalf_tmFS_R_dx_fn)
		  );
