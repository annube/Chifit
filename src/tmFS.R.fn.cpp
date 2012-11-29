

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


#include "utils.h"
#include "tmFS.R.fn.h"


#include <Rmath.h>




/*********************
 *
 * Evaluation of the R_0^k functions
 *
 *
 *********************/

typedef struct R_Integrand_Params_ {
  int k;
  double x,r;


} R_Integrand_Params;


complex<double> R_sigma(complex<double> z){
  return sqrt(1.-4./z);
}

complex<double> R_g(complex<double> z){
  return R_sigma(z) * log( ( R_sigma(z) - 1. ) / ( R_sigma(z) + 1. ) ) + 2.;
}



double R_integrand(double y, void *params){
  //  static int fnEvalCount=0;
  R_Integrand_Params *sp=(R_Integrand_Params *)params;
  //  ++fnEvalCount;
  //  cout << fnEvalCount << endl;


  if( sp->k % 2 == 0 ){
    return 
      pow(y,sp->k)
      * exp( -sp->x * sqrt( 1. + y * y )  )
      *
      R_g( complex<double>(2.,2.*sp->r*y) ).real();
  } else {
    return 
      pow(y,sp->k)
      * exp( -sp->x * sqrt( 1. + y * y )  )
      *
      R_g( complex<double>(2.,2.*sp->r*y) ).imag();
  }
    
}



double tmFS_R_fn(double k,double x,double r){

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
  F.function = &R_integrand;
  F.params = &rip;


  double intResult,intAbsErr;

  gsl_integration_qagi(&F,0,1.e-3,WS_SIZE,w,&intResult,&intAbsErr);

//   cout << "Integration result : " << intResult << endl;
//   cout << "Abs Error          : " << intAbsErr << endl;
//   cout << " intervals         : " << w->size << endl;

//  gsl_integration_workspace_free(w);

  return (intResult);

}



#include <Rcpp.h>

using namespace Rcpp;

RcppExport SEXP tmFS_R_fn_R(SEXP x,SEXP k, SEXP r){
  return wrap( tmFS_R_fn(as<int>(k),as<double>(x), as<double>(r) ) );
}
