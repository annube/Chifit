
#include <ginac/ginac.h>

#include <iostream>
#include <string>
#include <Rcpp.h>

#include <vector>
#include <map>

using namespace std;

using namespace GiNaC;


/// general method for evaluating an expression for a given set of parameters (parValues) and x values (xs)


typedef map<ex,vector<double> >::const_iterator XValueMapIt;

void eval_ex(ex Xpress, const exmap &parValues,const map<ex,vector<double> > &xs,vector<double> &result) {


    /* apply the map to the expression */
    ex Xpress_par_eval = evalf( Xpress.subs( parValues ) );



    int nx = ( xs.begin()->second ).size();

    for( int i =  0 ; i < nx ; i++){
      exmap xSubstMap;

      for( XValueMapIt it=xs.begin();it != xs.end() ; it++)
	xSubstMap[ it->first ] = ( it->second ) [ i ];

      ex_to<numeric>( 
		     Xpress_par_eval.subs( xSubstMap ) 
		     ) .to_double();

    }


  }
