

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



double R_integrand(double x, void *params){
  //  static int fnEvalCount=0;
  R_Integrand_Params *sp=(R_Integrand_Params *)params;
  //  ++fnEvalCount;
  //  cout << fnEvalCount << endl;
  return pow(x,sp->k)*exp(-sp->x*sqrt((1.+x*x/(sp->r*sp->r))))
    *
    R_g( complex<double>(2.,2.*x) ).real();

}



double tmFS_R_fn(double k,double x,double r){

//   cout << "k = " << (k) << endl
//        << "x = " << (x) << endl
//        << "r = " << (r) << endl ;

  int WS_SIZE=1000;

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

  gsl_integration_qagi(&F,1.e-14,1.e-3,WS_SIZE,w,&intResult,&intAbsErr);

//   cout << "Integration result : " << intResult << endl;
//   cout << "Abs Error          : " << intAbsErr << endl;
//   cout << " intervals         : " << w->size << endl;

//  gsl_integration_workspace_free(w);

  return (intResult);

}




