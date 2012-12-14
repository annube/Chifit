

#include <ginac/ginac.h>

using namespace GiNaC;


#include "utils.h"
#include "tmFS.dR.fn.h"
#include "tmFS.R.fn.h"
#include "tmFS.R.g.fn.h"




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


REGISTER_FUNCTION(tmFS_dR,
		  eval_func(eval_tmFS_dR_fn).
		  evalf_func(evalf_tmFS_dR_fn)
		  );
