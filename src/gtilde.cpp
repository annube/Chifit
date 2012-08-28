



#include <ginac/ginac.h>

using namespace GiNaC;

#include <iostream>
#include <iomanip>

using namespace std;

#include "Rmath.h"


#include "gtilde.h"
#include "gtilde.partitions.h"


const int GTILDE1_MAXORDER = 40;


double eval_gtilde1(double x){
  double sum=0.0;

  //  cout << x << endl;
  for( int i = 1 ; i <= GTILDE1_MAXORDER ; i ++ ){
    double sx = sqrt( (double) i ) * x;
    double ds = 4. * (double) getPartition(i) / ( sx ) * bessel_k( sx , 1. , 1.);
    sum += ds ; 
    //    cout << abs(ds/sum) << endl;
  }

  return sum;
}


double eval_dgtilde1(double x){
  double sum=0.0;

  //  cout << x << endl;
  for( int i = 1 ; i <= GTILDE1_MAXORDER ; i ++ ){
    double sx = sqrt( (double) i ) * x;
    double ds =
      4. * (double) getPartition(i) / ( sx ) *(
					       - 1. / x *  bessel_k( sx , 1. , 1.) 
					       -0.5 * sqrt( (double)i ) * ( bessel_k( sx , 0. , 1.)  + bessel_k( sx , 2. , 1.) ) );
;
    sum += ds ; 
    //    cout << abs(ds/sum) << endl;
  }

  return sum;
}


#include <Rcpp.h>


RcppExport SEXP eval_gtilde1_R( SEXP x){

  Rcpp::NumericVector vx(x);

  vx[0] = eval_gtilde1(vx[0]);

  return vx;
}


RcppExport SEXP eval_dgtilde1_R( SEXP x){

  Rcpp::NumericVector vx(x);

  vx[0] = eval_dgtilde1(vx[0]);

  return vx;
}

ex g_tilde_1_deri_eval(const ex &x){
  return dgtilde1(x).hold();
}


 ex g_tilde_1_deri_evalf(const ex &x){
  if( is_a<numeric>(x) ) {
    return eval_dgtilde1(ex_to<numeric>(x).to_double()) ;
  } else {
     return dgtilde1(x).hold();
  }
}


 ex g_tilde_1_deri_eval_deri(const ex &x,unsigned diffpar ){
  return 1./0.;
}


 void g_tilde_1_deri_print(const ex & arg, const print_context & c)
     {
         c.s << "gtilde_deri("; arg.print(c); c.s << ")";
     }





 ex g_tilde_1_eval(const ex &x){
  return gtilde1(x).hold();
}


 ex g_tilde_1_evalf(const ex &x){
  if( is_a<numeric>(x) ) {
    return eval_gtilde1(ex_to<numeric>(x).to_double()) ;
  } else {
     return gtilde1(x).hold();
  }
}


 ex g_tilde_1_eval_deri(const ex &x,unsigned diffpar ){
   //  cout << "called g_tilde_1_eval_deri" << endl;
     return dgtilde1(x);
 }



REGISTER_FUNCTION(dgtilde1, 
		  eval_func(g_tilde_1_deri_eval).
		  evalf_func(g_tilde_1_deri_evalf)//.
		  //		  derivative_func(g_tilde_1_deri_eval_deri).
// 		  latex_name("\\tilde{g}_1'").
// 		  print_func<print_dflt>(g_tilde_1_deri_print)
		  );

REGISTER_FUNCTION(gtilde1, 
		  eval_func(g_tilde_1_eval).
		  evalf_func(g_tilde_1_evalf).
 		  derivative_func(g_tilde_1_eval_deri));



