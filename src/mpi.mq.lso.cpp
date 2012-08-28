#include <ginac/ginac.h>

#include <iostream>
#include <string>
#include <Rcpp.h>

#include <vector>

using namespace std;

using namespace GiNaC;


#include "gtilde.h"
#include "eval.ex.lso.h"

namespace mpi_mq_lso{

  symbol aB0("aB_0");
  symbol af("af");
  symbol amu_q("a\\mu_q");
  symbol aLambda3("a\\Lambda_3");
  symbol aLambda4("a\\Lambda_4");
  symbol L("L");
  symbol ZP("Z_P");

  ex chi_mu_asq = 2.*aB0*amu_q/ZP; 
  ex xi_ll = chi_mu_asq / pow( 4. * Pi * af , 2. );



  typedef map< string  , vector<double> > AargsType;




  /**************************
   *
   *
   *        PION MASS
   *
   *
   *
   **************************/




  static ex getampisqXpression(){
    static ex X =  chi_mu_asq * ( 1.
					   + xi_ll  * log( chi_mu_asq / pow( aLambda3,2.) ) 
					   ) * pow(1. + 0.5 * xi_ll * gtilde1( sqrt(
										 chi_mu_asq ) *  L
									       ) ,2) ;
    return X;
  }



  RcppExport SEXP mpi_mq_lso(SEXP x, SEXP par,SEXP aargs,SEXP deri) {

    static  ex ampisq = getampisqXpression();

    Rcpp::NumericVector vpar(par);
    Rcpp::NumericVector vx(x);
    int nx = vx.size();

    Rcpp::NumericVector mpires(nx);
    Rcpp::NumericMatrix dmpires(nx,vpar.size());

    bool calcDeri = Rcpp::as<bool>(deri);

    Rcpp::List aargsMap ( aargs);
    int numLs = Rcpp::as<int>(aargsMap["numLs"]);
    Rcpp::NumericVector latSpacIndex = aargsMap["lsIndex"] ;
    Rcpp::NumericVector Loa = aargsMap["L"] ;
    Rcpp::NumericVector zpVal = aargsMap["ZP"] ;



    /// symbol for the ratio of two lattice spacings
    symbol R("R");
    int num_pure_parameters = vpar.size() - (numLs-1);

    vector<GiNaC::symbol> deriMap(4);
    deriMap[0] = aB0;
    deriMap[1] = af;
    deriMap[2] = aLambda3;
    deriMap[3] = R;


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



      ex X_R_subs = ampisq;

      /* 1.)  substitute R * par if not first lat. spac. */
      if( latSpacIndex[ix] > 1 ){
	exmap R_times_par;
	R_times_par[aB0] = R * aB0;
	R_times_par[af] = R * af;
	R_times_par[aLambda3] = R * aLambda3;
	X_R_subs = ampisq.subs( R_times_par );
      }


      /**
       * prepare numeric evaluation because it
       * can be used also for the deri 
       */
      exmap par_numeric_vals;
      /* parameters */
      par_numeric_vals[aB0] = vpar[0];
      par_numeric_vals[af] = vpar[1];
      par_numeric_vals[aLambda3] = vpar[2];
      par_numeric_vals[R] = vpar[3+latSpacIndex[ix]-2];

      /* regressors */
      par_numeric_vals[amu_q] = vx[ix];
      par_numeric_vals[ZP] = zpVal[ix];
      par_numeric_vals[L] = Loa[ix];



      if( calcDeri ){
      /* at this point we would do the derivatives */
	// for i_par in 1 to number of parameters
	int ipar = 0;
 	for( SymbolVecIt dit = deriMap.begin() ; dit != deriMap.end() ; dit++,ipar++){

	  if( (*dit) == R && latSpacIndex[ix] == 1 ) continue;

	  ex Xpress_deri = X_R_subs.diff( *dit ,1 );
	  ex Xpress_num_eval = Xpress_deri.subs( par_numeric_vals ).evalf();

	  if( is_a<numeric>(Xpress_num_eval) ){

	    if( (*dit) == R ) {
	      dmpires(ix,num_pure_parameters + latSpacIndex[ix]-2 ) = ex_to<numeric>( Xpress_num_eval ).to_double();
	    } else {
	      dmpires(ix,ipar) = ex_to<numeric>( Xpress_num_eval ).to_double();
	    }
	  } else {
	    cout << Xpress_num_eval << endl;
	  }
	}




      } else {

	/* perform numeric avaluation */
	ex Xpress_num_eval = X_R_subs.subs( par_numeric_vals ).evalf();
	
	if( is_a<numeric>(Xpress_num_eval) ){
	  mpires[ix] = ex_to<numeric>( Xpress_num_eval ).to_double();
	} else {
	  cout << Xpress_num_eval << endl;
	}
      }


    }


      if( calcDeri )
	return dmpires;
      else
	return mpires;

  }


  RcppExport SEXP mpi_mq_lso_gen(SEXP x, SEXP par,SEXP aargs,SEXP deri) {
    static  ex ampisq = getampisqXpression();

    /* the main parameters to optimize for */
    SymbolVec pureParVec;
    pureParVec.push_back(aB0);
    pureParVec.push_back(af);
    pureParVec.push_back(aLambda3);


    /* additional regressors besides the main regressor */
    SymbolStringVec ssvec;
    ssvec.push_back( SymbolStringPair( &L , "L" ) );
    ssvec.push_back( SymbolStringPair( &ZP , "ZP" ) );

    return    eval_ex_lso(ampisq, /* the expression to work on */
			  x,par,aargs,deri, /* pass on the parameters from R environment */
			  amu_q, /* the main regressor appearing in the expression */
			  pureParVec,  /* a vector of parameters to optimize for */
			  ssvec   /* a vector of additional regresssor and their name in the aargs list */
			  );
  }



};
