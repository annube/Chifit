
#include <ginac/ginac.h>

#include <iostream>
#include <string>
#include <Rcpp.h>

#include <vector>
#include <map>

using namespace std;

using namespace GiNaC;


#include "eval.ex.h"

/// general method for evaluating an expression for a given set of parameters (parValues) and x values (xs)
void eval_ex(ex Xpress, const exmap &parValues,const XValueMap &xs,Rcpp::NumericVector &result) {


  /* apply the map to the expression */
  ex Xpress_par_eval = evalf( Xpress.subs( parValues ) );


  int nx = ( xs.begin()->second )->size();

  for( int i =  0 ; i < nx ; i++){
    exmap xSubstMap;

    for( XValueMapIt it=xs.begin();it != xs.end() ; it++) {
      xSubstMap[ it->first ] =  (*( it->second )) [ i ];
    }

    result[i] = ex_to<numeric>( 
			       Xpress_par_eval.subs( xSubstMap )
			       ) .to_double();

  }


}


/// general method for evaluating the derivative of an expression for a given set of parameters (parValues) and x values (xs)

typedef exmap::const_iterator exmapCIt;
void eval_ex_deri(ex Xpress, const exmap &parValues,parList &smap,const XValueMap &xs,Rcpp::NumericMatrix &result) {

  int nx = ( xs.begin()->second )->size();
  
  int cpit=0;
  for( parListIt pit=smap.begin() ; pit != smap.end() ; pit++){
    
    ex dXpress_dpar = Xpress.diff( *(*pit) ,1);
    
    /* apply the map to the expression */
    ex dXpress_dpar_eval = evalf( dXpress_dpar.subs( parValues ) );
    
   
    for( int i =  0 ; i < nx ; i++){
      
      
      exmap xSubstMap;
      
      for( XValueMapIt it=xs.begin();it != xs.end() ; it++) {
	xSubstMap[ it->first ] =  (*( it->second )) [ i ];
	
      }
      
      result(i,cpit) = ex_to<numeric>( 
				      dXpress_dpar_eval.subs( xSubstMap )
				      ) .to_double();
      
    }
    ++cpit;
  }


}
