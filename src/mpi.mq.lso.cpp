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
  symbol Dmps("D_{m_{ps}}");
  symbol Dfps("D_{f_{ps}}");

  ex chi_mu_asq = 2. *  aB0 * amu_q/ZP; 
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
					   + Dmps) * pow(1. + 0.5 * xi_ll * gtilde1( sqrt(
										 chi_mu_asq ) *  L
									       ) ,2) ;
    return X;
  }


  RcppExport SEXP mpi_mq_lso(SEXP x, SEXP par,SEXP aargs,SEXP deri) {
    static  ex ampisq = getampisqXpression();

    /* the main parameters to optimize for */
    SymbolVec pureParVec;
    pureParVec.push_back(aB0);
    pureParVec.push_back(af);
    pureParVec.push_back(aLambda3);
    pureParVec.push_back(Dmps);

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
			  amu_q, /* the main regressor appearing in the expression */
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
    static ex X = af * ( 1 
			    - 2 *  xi_ll  * log( chi_mu_asq / pow( aLambda4,2) ) 
			 + Dfps
			    )
      *
      ( 1 - 2  * xi_ll * gtilde1( sqrt(  chi_mu_asq ) *  L    ) ) ;
    return X;
  }


  RcppExport SEXP fpi_mq_lso(SEXP x, SEXP par,SEXP aargs,SEXP deri) {
    static  ex afpi = getafpiXpression();

    /* the main parameters to optimize for */
    SymbolVec pureParVec;
    pureParVec.push_back(aB0);
    pureParVec.push_back(af);
    pureParVec.push_back(aLambda4);
    pureParVec.push_back(Dfps);

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
			  amu_q, /* the main regressor appearing in the expression */
			  pureParVec,  /* a vector of parameters to optimize for */
			  pureParDimE, /* dimension of parameter in terms of energy */
			  ssvec   /* a vector of additional regresssor and their name in the aargs list */
			  );
  }



};
