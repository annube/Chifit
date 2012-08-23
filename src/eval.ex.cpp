
#include <ginac/ginac.h>

#include <iostream>
#include <string>
#include <Rcpp.h>

#include <vector>
#include <map>

using namespace std;

using namespace GiNaC;


/// general method for evaluating an expression for a given set of parameters (parValues) and x values (xs)

void eval_ex(ex Xpress, const exmap &parValues,const map<ex,vector<double> > &xs) {


    /* apply the map to the expression */
    ex Xpress_par_eval = evalf( Xpress.subs( parValues ) );



    int nx = ( xs.begin()->second ).size();

    for( int i =  0 ; i < nx ; i++){
      exmap xSubstMap;
      xSubstMap[ xs.begin()->first ] = ( xs.begin()->second ) [ i ];

      // -->> now evaluat Xpress_par_eval numerically with xSubstMap

    }
    //
    //



//     for(int i = 0; i< vx.size() ; i++){
//       vres[i] = ex_to<numeric>( 
// 			       ampisq_num_eval.subs( lst( amu_q == vx[i] ,
// 							  a == latSpac[i] )) 
// 			       ) .
// 	to_double();
//     }


//     return vres;


  }
