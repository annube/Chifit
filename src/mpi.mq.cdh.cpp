







#include <ginac/ginac.h>

#include <iostream>
#include <string>
#include <Rcpp.h>

#include <vector>

using namespace std;

using namespace GiNaC;


#include "cdh.ex.h"
#include "eval.ex.lso.h"

namespace mpi_mq_cdh{

  symbol B0("B_0");
  symbol f0("f");
  symbol mu_q("\\mu_q");
  symbol Lambda1("\\Lambda_1");
  symbol Lambda2("\\Lambda_2");
  symbol Lambda3("\\Lambda_3");
  symbol Lambda4("\\Lambda_4");
  symbol L("L");
  symbol ZP("Z_P");
  symbol Dmps("D_{m_{ps}}");
  symbol Dfps("D_{f_{ps}}");

  ex chi_mu_sq = 2. * B0 * mu_q/ZP; 
  ex xi_ll = chi_mu_sq / pow( 4. * Pi * f0 , 2. );



  typedef map< string  , vector<double> > AargsType;



  /**************************
   *
   *
   *        PION MASS
   *
   *
   *
   **************************/


  static ex getmpisqXpression(){
    static ex X =  chi_mu_sq * ( 1.
					   + xi_ll  * log( chi_mu_sq / pow( Lambda3,2.) ) 
				 ) 
      * 
      pow( 
	  cdh_mpi_ex( chi_mu_sq , f0 , L , Lambda1,Lambda2,Lambda3,Lambda4 )
	  , 2
	  ) ;
    return X;
  }


  RcppExport SEXP mpi_mq_cdh(SEXP x, SEXP par,SEXP aargs,SEXP deri) {
    static  ex mpisq = getmpisqXpression();

    /* the main parameters to optimize for */
    SymbolVec pureParVec;
    pureParVec.push_back(B0);
    pureParVec.push_back(f0);
    pureParVec.push_back(Lambda1);
    pureParVec.push_back(Lambda2);
    pureParVec.push_back(Lambda3);
    pureParVec.push_back(Lambda4);

    vector<double> pureParDimE;
    pureParDimE.push_back(1.);
    pureParDimE.push_back(1.);
    pureParDimE.push_back(1.);
    pureParDimE.push_back(1.);
    pureParDimE.push_back(1.);
    pureParDimE.push_back(1.);


    /* additional regressors besides the main regressor */
    SymbolStringVec ssvec;
    ssvec.push_back( SymbolStringPair( &L , "L" ) );
    ssvec.push_back( SymbolStringPair( &ZP , "ZP" ) );

    return    eval_ex_lso(mpisq, /* the expression to work on */
			  x,par,aargs,deri, /* pass on the parameters from R environment */
			  mu_q, /* the main regressor appearing in the expression */
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


  static ex getfpiXpression(){
    static ex X = f0 * ( 1 
			 - 2. *  xi_ll  * log( chi_mu_sq / pow( Lambda4,2) ) 
			 )
      *
      cdh_fpi_ex( chi_mu_sq , f0 , L , Lambda1,Lambda2,Lambda3,Lambda4 ) ;
    return X;
  }



  RcppExport SEXP fpi_mq_cdh(SEXP x, SEXP par,SEXP aargs,SEXP deri) {
    static  ex mpisq = getfpiXpression();

    /* the main parameters to optimize for */
    SymbolVec pureParVec;
    pureParVec.push_back(B0);
    pureParVec.push_back(f0);
    pureParVec.push_back(Lambda1);
    pureParVec.push_back(Lambda2);
    pureParVec.push_back(Lambda3);
    pureParVec.push_back(Lambda4);

    vector<double> pureParDimE;
    pureParDimE.push_back(1.);
    pureParDimE.push_back(1.);
    pureParDimE.push_back(1.);
    pureParDimE.push_back(1.);
    pureParDimE.push_back(1.);
    pureParDimE.push_back(1.);


    /* additional regressors besides the main regressor */
    SymbolStringVec ssvec;
    ssvec.push_back( SymbolStringPair( &L , "L" ) );
    ssvec.push_back( SymbolStringPair( &ZP , "ZP" ) );

    return    eval_ex_lso(mpisq, /* the expression to work on */
			  x,par,aargs,deri, /* pass on the parameters from R environment */
			  mu_q, /* the main regressor appearing in the expression */
			  pureParVec,  /* a vector of parameters to optimize for */
			  pureParDimE,
			  ssvec   /* a vector of additional regresssor and their name in the aargs list */
			  );
  }


};
