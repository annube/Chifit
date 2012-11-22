
#include <math.h>

#include <ginac/ginac.h>

using namespace GiNaC;

#include <utility>
#include <iostream>

using namespace std;


#include "eval.ex.lso.h"


void debug_print_symbol( symbol &s){ cout << s << endl;}

double numEvalXPression(ex Xpression,exmap substitutions){
    /* perform numeric avaluation */
    ex Xpress_num_eval = Xpression.subs( substitutions ).evalf();
    
    
    if( is_a<numeric>(Xpress_num_eval) ){
      return ex_to<numeric>( Xpress_num_eval ).to_double();
    } else {
      cout << "ERROR: could not evaluate expression numerically."<< endl
	   << "This is the remaining expression (most likly with unresolved symbols): " << endl;
      cout << Xpress_num_eval << endl;
      return NAN;
    }
}


enum ParameterGroup {PG_Scaling_Par=0,PG_R,PG_LS_Dep};

/**
 * for definiteness we have to fix the order of the parameters
 *
 * 1.) Parameters scaling with the lattice spacing
 * 2.) the ratios of the lattice spacings
 * 3.) lattice spacing dependent (but not in a trivial way, not scaling with R) parameters
 *
 */

SEXP eval_ex_lso(
		 ex Xpression,SEXP x,SEXP par, SEXP aargs,SEXP deri,
		 symbol mainRegressor,
		 SymbolVec pureParVec,
		 vector<double> pureParDimE,
		 SymbolStringVec addRegr,
		 SymbolVec lsDepPar
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


  int num_pure_parameters = pureParVec.size();

  /* deri map will contain all symbols that the expression will be derived for */


  /* vector containing the kind of group to which the parameters belongs*/
  vector<ParameterGroup> deriParGroup;
  
  /* 1.) the pure parameters  */
  SymbolVec deriMap(pureParVec);
  // -> which parameter group  ?!
  deriParGroup.insert( deriParGroup.end(),pureParVec.size(),PG_Scaling_Par);

  /* 2.) the ratio to the largest lattice spacing */
  deriMap.push_back(R);
  // -> which parameter group  ?!
  deriParGroup.insert( deriParGroup.end(),1,PG_R);

  /* 3.) parameters depending on the lattice spacing but do not scale trivially with the lattice spacing */
  deriMap.insert(deriMap.end(),lsDepPar.begin(),lsDepPar.end());
  // -> which parameter group  ?!
  deriParGroup.insert( deriParGroup.end(),lsDepPar.size(),PG_LS_Dep);


  ex X_R_subs = Xpression;
  /* 1.)  substitute R * par if not first lat. spac. */
  exmap R_times_par;
  
  for( int pi = 0 ; pi < pureParVec.size() ; pi++)
    R_times_par[ pureParVec[pi] ]  = pureParVec[pi] * pow( R , pureParDimE[pi] );
    
  X_R_subs = Xpression.subs( R_times_par );


  ExVec Xpressions_to_evaluate;


  if( ! calcDeri){
    Xpressions_to_evaluate.push_back( X_R_subs );
  } else  {
    /* push back the derivatives of X_R_subs */
    for(unsigned long ipar = 0 ; ipar < deriMap.size() ; ipar++){
      Xpressions_to_evaluate.push_back( X_R_subs.diff( deriMap[ipar] ) );
    }
  }


  /**
   * loop over all observations 
   */
  for(int ix=0;ix<nx;ix++){

    /* first create the substitution map */
    exmap par_numeric_vals;


    /* we have to replace five kinds of symbols in the expression: */

    /* 1.) the ordinary parameters f0, B0, c2, ...*/
    int lin_count = 0;
    for( SymbolVecIt svit = pureParVec.begin() ; svit!=pureParVec.end() ; svit++,lin_count++) {
      par_numeric_vals[*svit] = vpar[lin_count];
    }


    /* 2.) the lattice spacing scaling parameter R */
    if( latSpacIndex[ix] > 1 ) {
      par_numeric_vals[R] = vpar[num_pure_parameters+latSpacIndex[ix]-2];
    } else {
      par_numeric_vals[R] = 1.;
    }


    /* 3.) regressors \mu_q (unrenormalized q-mass)*/
    par_numeric_vals[mainRegressor] = vx[ix];

    /* 4.) additional regressors e.g. L/a */
    for( SymbolVecVecIt svvIt = addRegressorsValues.begin() ; svvIt!=addRegressorsValues.end() ; svvIt++) {
      par_numeric_vals[ *(svvIt->first)] = (svvIt->second)[ix];
    }


    /* 5.) lattice spacing dependent parameters (not via R) mainly Z_P*/
    for( SymbolVecIt lsDepParIt=lsDepPar.begin(); lsDepParIt != lsDepPar.end(); lsDepParIt++) {
      int index=num_pure_parameters + (numLs-1) + latSpacIndex[ix]-1;
      par_numeric_vals[ *lsDepParIt ] = vpar[index];
    }
  

    if( ! calcDeri ){  evalres[ix] = numEvalXPression(Xpressions_to_evaluate[0],par_numeric_vals); }
    else {
      /* loop over expressions */
      for( unsigned long  expI = 0 ; expI < Xpressions_to_evaluate.size() ; expI++){
	//derires(ix,num_pure_parameters + latSpacIndex[ix]-2 ) = ex_to<numeric>( Xpress_num_eval ).to_double();


	switch( deriParGroup[expI] ){

	case PG_Scaling_Par : 
	  derires(ix,expI) = numEvalXPression(Xpressions_to_evaluate[expI],par_numeric_vals);
	  break;

	case PG_R : 
	  if( latSpacIndex[ix] > 1 ){
	    derires(ix,num_pure_parameters+latSpacIndex[ix]-2) = numEvalXPression(Xpressions_to_evaluate[expI],par_numeric_vals);
	  }
	  break;

	case PG_LS_Dep :
	  derires(ix,num_pure_parameters+numLs-1+latSpacIndex[ix]-1) = numEvalXPression(Xpressions_to_evaluate[expI],par_numeric_vals);
	  break;
	default: cerr << "Error: We have a parameter belonging to an unkown parameter group!!" << endl; break;
	}

      }
    }

  }


  if( ! calcDeri )
    return( evalres );
  else
    return (derires);


}
