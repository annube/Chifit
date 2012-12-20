

#include <ginac/ginac.h>

using namespace GiNaC;


#include <Rcpp.h>

using namespace Rcpp;

#include "utils.h"
#include "tmFS.dR.fn.h"
#include "tmFS.R.fn.h"
#include "tmFS.R.g.fn.h"



/**
 * integrand for the R function with g' instead of g
 */
double dR_integrand(double y,void *params){
  //  static int fnEvalCount=0;
  R_Integrand_Params *sp=(R_Integrand_Params *)params;
  //  ++fnEvalCount;
  //  cout << fnEvalCount << endl;


  if( sp->k % 2 == 0 ){
    return 
      pow(y,sp->k)
      * exp( -sp->x * sqrt( sp->n* ( 1. + y * y ))  )
      *
      chifit::R_dg( complex<double>(2.,2.*sp->r*y) ).real();
  } else {
    return 
      pow(y,sp->k)
      * exp( -sp->x * sqrt( sp->n* (1. + y * y) )  )
      *
      chifit::R_dg( complex<double>(2.,2.*sp->r*y) ).imag();
  }
}

double dR_integrand_dx(double y, void *params){
  //  static int fnEvalCount=0;
  R_Integrand_Params *sp=(R_Integrand_Params *)params;
  //  ++fnEvalCount;
  //  cout << fnEvalCount << endl;


  if( sp->k % 2 == 0 ){
    return 
      - pow(y,sp->k) * sqrt( sp->n*(1. + y * y) )
      * exp( -sp->x * sqrt( sp->n*(1. + y * y ))  )
      *
      chifit::R_dg( complex<double>(2.,2.*sp->r*y) ).real();
  } else {
    return 
      - pow(y,sp->k) * sqrt( sp->n* ( 1. + y * y ) )
      * exp( -sp->x * sqrt( sp->n* ( 1. + y * y )  ) )
      *
      chifit::R_dg( complex<double>(2.,2.*sp->r*y) ).imag();
  }
    
}


double dR_integrand_dr(double y, void *params){
  //  static int fnEvalCount=0;
  R_Integrand_Params *sp=(R_Integrand_Params *)params;
  //  ++fnEvalCount;
  //  cout << fnEvalCount << endl;


  if( ( sp->k ) % 2 == 0 ){
    return 
      pow(y,sp->k) 
      * exp( -sp->x * sqrt( sp->n* ( 1. + y * y )  ) ) 
      *
      ( -2. ) * y *( chifit::R_ddg( complex<double>(2.,2.*sp->r*y) ) ).imag();
  } else {
    return 
      pow(y,sp->k)
      * exp( -sp->x * sqrt( sp->n* ( 1. + y * y ) ) )
      *
      ( chifit::R_ddg( complex<double>(2.,2.*sp->r*y) ) ).real()* 2. * y;
  }
    
}




/**
 * ginac function for evaluation
 */

ex eval_tmFS_dR_fn(const ex &x,const ex & k , const ex & r){
  if( is_a<numeric>(k) ) 
    return tmFS_dR(x,ex_to<numeric>(k).to_int(),r).hold();
  else 
    return tmFS_dR(x,k,r).hold();
}


ex evalf_tmFS_dR_fn(const ex &x,const ex & k , const ex & r){
  if( is_a<numeric>(x) &&  is_a<numeric>(k) && is_a<numeric>(r) ) {

    return eval_tmFS_R(
		      ex_to<numeric>(x).to_double(),
		      MyExToInt( k ),
		      ex_to<numeric>(r).to_double(),
		      &dR_integrand
		      ) ;
  } else {
    return tmFS_R(x,k,r).hold();
  }
}

ex eval_tmFS_dR_fn_deri(const ex &x,const ex & k , const ex & r,unsigned diffpar){
  if(diffpar == 0 )
    return tmFS_dR_dx(x,k,r);
  else if( diffpar == 2 ){
    return tmFS_dR_dr(x,k,r);
  } else {
    cerr << "Error : dR can only be derived w.r.t. x and r " << endl;
  }
  return NAN;
}


/**
 * make this function be known to ginac
 */

REGISTER_FUNCTION(tmFS_dR,
		  eval_func(eval_tmFS_dR_fn).
		  evalf_func(evalf_tmFS_dR_fn).
		  derivative_func(eval_tmFS_dR_fn_deri)
		  );



/**
 * export the function to R
 */
RcppExport SEXP tmFS_dR_fn_R_ginac(SEXP x,SEXP k, SEXP r){
  ex fnv=tmFS_dR(as<double>(x),as<int>(k),as<double>(r));
  return wrap( ex_to<numeric>(fnv.evalf()).to_double() );
}



/**
 * ginac function for evaluation
 */

ex eval_tmFS_dR_dx_fn(const ex &x,const ex & k , const ex & r){
  if( is_a<numeric>(k) ) 
    return tmFS_dR_dx(x,ex_to<numeric>(k).to_int(),r).hold();
  else 
    return tmFS_dR_dx(x,k,r).hold();
}


ex evalf_tmFS_dR_dx_fn(const ex &x,const ex & k , const ex & r){
  if( is_a<numeric>(x) &&  is_a<numeric>(k) && is_a<numeric>(r) ) {

    return eval_tmFS_R(
		      ex_to<numeric>(x).to_double(),
		      MyExToInt( k ),
		      ex_to<numeric>(r).to_double(),
		      &dR_integrand_dx
		      ) ;
  } else {
    return tmFS_R_dx(x,k,r).hold();
  }
}


REGISTER_FUNCTION(tmFS_dR_dx,
		  eval_func(eval_tmFS_dR_dx_fn).
		  evalf_func(evalf_tmFS_dR_dx_fn)
		  );



RcppExport SEXP tmFS_dR_dx_fn_R_ginac(SEXP x,SEXP k, SEXP r){
  ex fnv=tmFS_dR_dx(as<double>(x),as<int>(k),as<double>(r));
  return wrap( ex_to<numeric>(fnv.evalf()).to_double() );
}

/**
 * ginac function for evaluation
 */

ex eval_tmFS_dR_dr_fn(const ex &x,const ex & k , const ex & r){
  if( is_a<numeric>(k) ) 
    return tmFS_dR_dr(x,ex_to<numeric>(k).to_int(),r).hold();
  else 
    return tmFS_dR_dr(x,k,r).hold();
}


ex evalf_tmFS_dR_dr_fn(const ex &x,const ex & k , const ex & r){
  if( is_a<numeric>(x) &&  is_a<numeric>(k) && is_a<numeric>(r) ) {

    return eval_tmFS_R(
		      ex_to<numeric>(x).to_double(),
		      MyExToInt( k ),
		      ex_to<numeric>(r).to_double(),
		      &dR_integrand_dr
		      ) ;
  } else {
    return tmFS_dR_dr(x,k,r).hold();
  }
}


REGISTER_FUNCTION(tmFS_dR_dr,
		  eval_func(eval_tmFS_dR_dr_fn).
		  evalf_func(evalf_tmFS_dR_dr_fn)
		  );



RcppExport SEXP tmFS_dR_dr_fn_R_ginac(SEXP x,SEXP k, SEXP r){
  ex fnv=tmFS_dR_dr(as<double>(x),as<int>(k),as<double>(r));
  return wrap( ex_to<numeric>(fnv.evalf()).to_double() );
}
