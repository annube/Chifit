

#include <ginac/ginac.h>

using namespace GiNaC;

#include <utility>
#include <iostream>

using namespace std;


#include "eval.ex.lso.h"


void debug_print_symbol( symbol &s){ cout << s << endl;}


SEXP eval_ex_lso(
		 ex Xpression,SEXP x,SEXP par, SEXP aargs,SEXP deri,
		 symbol mainRegressor,
		 SymbolVec pureParVec,
		 vector<double> pureParDimE,
		 SymbolStringVec addRegr
		 ) {

  Rcpp::NumericVector vpar(par);
  Rcpp::NumericVector vx(x);
  int nx = vx.size();

  Rcpp::NumericVector evalres(nx);
  Rcpp::NumericMatrix derires(nx,vpar.size());

  bool calcDeri = Rcpp::as<bool>(deri);

  Rcpp::List aargsMap ( aargs);
  int numLs = Rcpp::as<int>(aargsMap["numLs"]);
  Rcpp::NumericVector latSpacIndex = aargsMap["lsIndex"] ;



  SymbolVecVec addRegressorsValues;


  for( SymbolStringVecIt ssvIt= addRegr.begin() ; ssvIt != addRegr.end(); ssvIt++)
    addRegressorsValues.push_back( SymbolVecPair( (ssvIt->first) , aargsMap[ssvIt->second] ) );




  /// symbol for the ratio of two lattice spacings
  symbol R("R");


  int num_pure_parameters = vpar.size() - (numLs-1);

  /* deri map will contain all symbols that the expression will be derived for */
  SymbolVec deriMap(pureParVec);
  deriMap.push_back(R);

  /* loop over all observations (dimension of one regressor vector) */
  for(int ix=0;ix<nx;ix++){

    /**
     *  pseudocode:
     *  1.) if we are not at the first lattice spacing
     *      replace the parameters by R * the Parametere
     *  2.) if we calculate the derivative derive Xpress.
     *      with respect to the parameter under consideration
     *  3.) substitute the numerical values of the parameter,
     *      the quark mass and other values like L and ZP
     */



    ex X_R_subs = Xpression;

    /* 1.)  substitute R * par if not first lat. spac. */
    if( latSpacIndex[ix] > 1 ){

      exmap R_times_par;


      for( int pi = 0 ; pi < pureParVec.size() ; pi++)
	R_times_par[ pureParVec[pi] ]  = pureParVec[pi] * pow( R , pureParDimE[pi] );


      X_R_subs = Xpression.subs( R_times_par );
    }



    /**
     * prepare numeric evaluation because it
     * can be used also for the deri 
     */
    exmap par_numeric_vals;
    /* parameters */
    int lin_count = 0;
    for( SymbolVecIt svit = pureParVec.begin() ; svit!=pureParVec.end() ; svit++,lin_count++)
      par_numeric_vals[*svit] = vpar[lin_count];


    if( latSpacIndex[ix] > 1 ) 
      par_numeric_vals[R] = vpar[num_pure_parameters+latSpacIndex[ix]-2];


    /* regressors */
    par_numeric_vals[mainRegressor] = vx[ix];

    /* additional regressors */
    for( SymbolVecVecIt svvIt = addRegressorsValues.begin() ; svvIt!=addRegressorsValues.end() ; svvIt++)
      par_numeric_vals[ *(svvIt->first)] = (svvIt->second)[ix];


    if( calcDeri ){
      /* at this point we would do the derivatives */
      // for i_par in 1 to number of parameters
      int ipar = 0;
      for( SymbolVecIt dit = deriMap.begin() ; dit != deriMap.end() ; dit++,ipar++){

	if( (*dit) == R && latSpacIndex[ix] == 1 ) continue;

	ex Xpress_deri = X_R_subs.diff( (*dit) ,1 );

	ex Xpress_num_eval = Xpress_deri.subs( par_numeric_vals ).evalf();

	if( is_a<numeric>(Xpress_num_eval) ){

	  if( (*dit) == R ) {
	    derires(ix,num_pure_parameters + latSpacIndex[ix]-2 ) = ex_to<numeric>( Xpress_num_eval ).to_double();
	  } else {
	    derires(ix,ipar) = ex_to<numeric>( Xpress_num_eval ).to_double();
	  }
	} else {
	  cout << Xpress_num_eval << endl;
	}
      }




    } else {

      /* perform numeric avaluation */
      ex Xpress_num_eval = X_R_subs.subs( par_numeric_vals ).evalf();
	

      if( is_a<numeric>(Xpress_num_eval) ){
	evalres[ix] = ex_to<numeric>( Xpress_num_eval ).to_double();
      } else {
	cout << Xpress_num_eval << endl;
      }
    }


  }


  if( calcDeri )
    return derires;
  else
    return evalres;

}
