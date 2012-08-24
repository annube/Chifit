
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


    //    cout << Xpress_par_eval << endl;


     int nx = ( xs.begin()->second )->size();

    for( int i =  0 ; i < nx ; i++){
      exmap xSubstMap;

      for( XValueMapIt it=xs.begin();it != xs.end() ; it++) {
	xSubstMap[ it->first ] =  (*( it->second )) [ i ];

	//	cout << "substituting " << it->first << " = " << (*( it->second )) [ i ] << endl;
      }

      //      cout << Xpress_par_eval.subs( xSubstMap ) << endl;

      result[i] = ex_to<numeric>( 
		     Xpress_par_eval.subs( xSubstMap )
		     ) .to_double();

    }


  }


/// general method for evaluating an expression for a given set of parameters (parValues) and x values (xs)

typedef exmap::const_iterator exmapCIt;
void eval_ex_deri(ex Xpress, const exmap &parValues,parList &smap,const XValueMap &xs,Rcpp::NumericMatrix &result) {




    //    cout << Xpress_par_eval << endl;


     int nx = ( xs.begin()->second )->size();

     cout << "nx = " << nx << endl;

     int cpit=0;
     for( parListIt pit=smap.begin() ; pit != smap.end() ; pit++){

       cout << "cpit=" << cpit << endl;
       ex dXpress_dpar = Xpress.diff( *(*pit) ,1);

       /* apply the map to the expression */
       ex dXpress_dpar_eval = evalf( dXpress_dpar.subs( parValues ) );

       cout << dXpress_dpar_eval << endl;

       for( int i =  0 ; i < nx ; i++){


	 exmap xSubstMap;
 
	 for( XValueMapIt it=xs.begin();it != xs.end() ; it++) {
	   xSubstMap[ it->first ] =  (*( it->second )) [ i ];
	   
	   //	cout << "substituting " << it->first << " = " << (*( it->second )) [ i ] << endl;
	 }
	 
	 //      cout << Xpress_par_eval.subs( xSubstMap ) << endl;
	 
	 result[i,cpit] = ex_to<numeric>( 
				      dXpress_dpar_eval.subs( xSubstMap )
				      ) .to_double();
	 
       }
       ++cpit;
     }


  }
