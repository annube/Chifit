#include <ginac/ginac.h>

#include <iostream>
#include <string>
#include <Rcpp.h>

#include <vector>

using namespace std;

using namespace GiNaC;


#include "gtilde.h"

#include "eval.ex.lso.h"
#include "symbols.h"

namespace chifit{


  ex chi_mu_asq = 2. *  B * mu/ZP; 



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
    ex xi_ll = chi_mu_asq / pow( 4. * Pi * f , 2. );
    static ex X =  chi_mu_asq * ( 1.
					   + xi_ll  * log( chi_mu_asq / pow( Lambda3,2.) ) 
					   + CMpm) * pow(1. + 0.5 * xi_ll * gtilde1( sqrt(
										 chi_mu_asq ) *  L
									       ) ,2.) ;
    return X;
  }


  RcppExport SEXP mpi_mq_lso(SEXP x, SEXP par,SEXP aargs,SEXP deri) {
    static  ex ampisq = getampisqXpression();

    /* the main parameters to optimize for */
    SymbolVec pureParVec;
    pureParVec.push_back(B);
    pureParVec.push_back(f);
    pureParVec.push_back(Lambda3);
    pureParVec.push_back(CMpm);

    vector<double> pureParDimE;
    pureParDimE.push_back(1.);
    pureParDimE.push_back(1.);
    pureParDimE.push_back(1.);
    pureParDimE.push_back(2.);


    /* additional regressors besides the main regressor */
    SymbolStringVec ssvec;
    ssvec.push_back( SymbolStringPair( &L , "L" ) );
    ssvec.push_back( SymbolStringPair( &ZP , "ZP" ) );

    return    eval_ex_lso(ampisq, /* the expression to work on */
			  x,par,aargs,deri, /* pass on the parameters from R environment */
			  mu, /* the main regressor appearing in the expression */
			  pureParVec,  /* a vector of parameters to optimize for */
			  pureParDimE,
			  ssvec   /* a vector of additional regresssor and their name in the aargs list */
			  );
  }


  /**************************
   *
   *
   *   PION DECAY CONSTANT
   *
   *
   *
   **************************/


  static ex getafpiXpression(){
    ex xi_ll = chi_mu_asq / pow( 4. * Pi * f , 2. );
    static ex X = f * ( 1 
			    - 2. *  xi_ll  * log( chi_mu_asq / pow( Lambda4,2) ) 
			 + Cf
			    )
      *
      ( 1. - 2.  * xi_ll * gtilde1( sqrt(  chi_mu_asq ) *  L    ) ) ;
    return X;
  }


  RcppExport SEXP fpi_mq_lso(SEXP x, SEXP par,SEXP aargs,SEXP deri) {
    static  ex afpi = getafpiXpression();

    /* the main parameters to optimize for */
    SymbolVec pureParVec;
    pureParVec.push_back(B);
    pureParVec.push_back(f);
    pureParVec.push_back(Lambda4);
    pureParVec.push_back(Cf);

    vector<double> pureParDimE;
    pureParDimE.push_back(1.);
    pureParDimE.push_back(1.);
    pureParDimE.push_back(1.);
    pureParDimE.push_back(2.);



    /* additional regressors besides the main regressor */
    SymbolStringVec ssvec;
    ssvec.push_back( SymbolStringPair( &L , "L" ) );
    ssvec.push_back( SymbolStringPair( &ZP , "ZP" ) );

    return    eval_ex_lso(afpi, /* the expression to work on */
			  x,par,aargs,deri, /* pass on the parameters from R environment */
			  mu, /* the main regressor appearing in the expression */
			  pureParVec,  /* a vector of parameters to optimize for */
			  pureParDimE, /* dimension of parameter in terms of energy */
			  ssvec   /* a vector of additional regresssor and their name in the aargs list */
			  );
  }



};
