#include <ginac/ginac.h>

#include <iostream>
#include <string>
#include <Rcpp.h>

#include <vector>

using namespace std;
using namespace GiNaC;


#include "gtilde.h"
#include "eval.ex.lso.h"

#include "mpi.mq.ob.h"


/* we need a namespace here because the same symbols as below are defined in mpi.mq.cpp as well */
namespace mpi_mq_ob {

  symbol B("aB_0");
  symbol f("af");
  symbol mu("a\\mu_q");
  symbol Lambda3("a\\Lambda_3");
  symbol Lambda4("a\\Lambda_4");
  symbol Xi3("\\Xi_3");
  symbol L("L");
  symbol ZP("Z_P");
  symbol CMpm("C_{M_{\\pm}}");
  symbol Cf("D_{f}");
  symbol CM0("C_{M_0}");
  symbol c2("c_2");

  ex Mpm_sq = 2*B*mu/ZP; 
  ex M0_sq = sqrt( pow(  2*B*mu/ZP + 2*c2 ,2 ) ); 


//   ex xi_pm = Mpm_sq/ pow( 4. * Pi * f , 2 );
//   ex xi_0  = M0_sq/ pow( 4. * Pi * f , 2 );


//   ex fse_log_corr_Mpm = gtilde1( sqrt( Mpm_sq ) *  L );
//   ex fse_log_corr_M0  = gtilde1( sqrt( M0_sq ) *  L );

  
//   ex log_M0_L3 = log( M0_sq / pow( Lambda3 , 2 ) ) + fse_log_corr_M0 ;
//   ex log_M0_L4 = log( M0_sq / pow( Lambda4 , 2 ) ) + fse_log_corr_M0;

//   ex log_Mpm_L3 = log( Mpm_sq / pow( Lambda3 , 2 ) ) + fse_log_corr_Mpm;
//   ex log_Mpm_L4 = log( Mpm_sq / pow( Lambda4 , 2 ) ) + fse_log_corr_Mpm;


  ex get_M_pm_sq_Xpression(){
    static ex xi_0  = M0_sq/ pow( 4. * Pi * f , 2 );
    static ex fse_log_corr_M0  = gtilde1( sqrt( M0_sq ) *  L );
    static ex log_M0_L3 = log( M0_sq / pow( Lambda3 , 2 )  ) + fse_log_corr_M0 ;

    static ex X = (  Mpm_sq * ( 1. + xi_0 * log_M0_L3 + CMpm ) ) ; 
    return X;
  }


  ex get_M_0_sq_Xpression(){

    static ex xi_0  = M0_sq/ pow( 4. * Pi * f , 2 );
    static ex fse_log_corr_M0  = gtilde1( sqrt( M0_sq ) *  L );
    static ex log_M0_L3 = log( M0_sq / pow( Lambda3 , 2 )  ) + fse_log_corr_M0 ;

    static ex xi_pm  = Mpm_sq/ pow( 4. * Pi * f , 2 );
    static ex fse_log_corr_Mpm  = gtilde1( sqrt( Mpm_sq ) *  L );
    static ex log_Mpm_L3 = log( Mpm_sq / pow( Lambda3 , 2 )  ) + fse_log_corr_Mpm ;

    static ex log_M0_X3 = log( M0_sq / pow( Xi3 , 2 )  ) + fse_log_corr_M0 ;

    static ex X =   Mpm_sq * ( 1. +2.*xi_pm * log_Mpm_L3 - xi_0 * log_M0_L3 ) 
      + 2 * c2 *( 1 - 4. * xi_0 * log_M0_X3 + CM0); 
    return X;
  }


  ex get_f_pm_Xpression(){
    static ex xi_pm = Mpm_sq/ pow( 4. * Pi * f , 2 );
    static ex xi_0  = M0_sq/ pow( 4. * Pi * f , 2 );
    static ex fse_log_corr_Mpm = gtilde1( sqrt( Mpm_sq ) *  L );
    static ex fse_log_corr_M0  = gtilde1( sqrt( M0_sq ) *  L );
    static ex log_M0_L4 = log( M0_sq / pow( Lambda4 , 2 )  ) + fse_log_corr_M0;
    static ex log_Mpm_L4 = log( Mpm_sq / pow( Lambda4 , 2 ) ) + fse_log_corr_Mpm;

    static ex X =  f * ( 1. - ( xi_pm * log_Mpm_L4 + xi_0 * log_M0_L4 ) + Cf ) ;
    return X;
  }




  /**************************
   *
   *
   *    CHARGED PION MASS
   *
   *
   *
   **************************/

  RcppExport SEXP mpi_mq_ob(SEXP x, SEXP par,SEXP aargs,SEXP deri) {
    static  ex mpisq = get_M_pm_sq_Xpression();

    /* the main parameters to optimize for */
    SymbolVec pureParVec;
    vector<double> pureParDimE;

    /* the parameters to be fitted and their energy dimension */

    pureParVec.push_back(B);       pureParDimE.push_back(1.);
    pureParVec.push_back(f);       pureParDimE.push_back(1.);
    pureParVec.push_back(c2);      pureParDimE.push_back(4.);
    pureParVec.push_back(Lambda3); pureParDimE.push_back(1.);
    pureParVec.push_back(CMpm);    pureParDimE.push_back(2.);

    


    /* additional regressors besides the main regressor */
    SymbolStringVec ssvec;
    ssvec.push_back( SymbolStringPair( &L , "L" ) );
    ssvec.push_back( SymbolStringPair( &ZP , "ZP" ) );

    return    eval_ex_lso(mpisq, /* the expression to work on */
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
   *    NEUTRAL PION MASS
   *
   *
   *
   **************************/

  RcppExport SEXP mpi_0_mq_ob(SEXP x, SEXP par,SEXP aargs,SEXP deri) {
    static  ex mpisq = get_M_0_sq_Xpression();

    /* the main parameters to optimize for */
    SymbolVec pureParVec;
    vector<double> pureParDimE;

    /* the parameters to be fitted and their energy dimension */

    pureParVec.push_back(B);       pureParDimE.push_back(1.);
    pureParVec.push_back(f);       pureParDimE.push_back(1.);
    pureParVec.push_back(c2);      pureParDimE.push_back(4.);
    pureParVec.push_back(Lambda3); pureParDimE.push_back(1.);
    pureParVec.push_back(Xi3);     pureParDimE.push_back(1.);
    pureParVec.push_back(CM0);    pureParDimE.push_back(2.);

    


    /* additional regressors besides the main regressor */
    SymbolStringVec ssvec;
    ssvec.push_back( SymbolStringPair( &L , "L" ) );
    ssvec.push_back( SymbolStringPair( &ZP , "ZP" ) );

    return    eval_ex_lso(mpisq, /* the expression to work on */
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
   *    PION DECAY
   *
   *
   *
   **************************/

  RcppExport SEXP fpi_mq_ob(SEXP x, SEXP par,SEXP aargs,SEXP deri) {
    static  ex fpisq = get_f_pm_Xpression();

    /* the main parameters to optimize for */
    SymbolVec pureParVec;
    vector<double> pureParDimE;

    /* the parameters to be fitted and their energy dimension */

    pureParVec.push_back(B);       pureParDimE.push_back(1.);
    pureParVec.push_back(f);       pureParDimE.push_back(1.);
    pureParVec.push_back(c2);      pureParDimE.push_back(4.);
    pureParVec.push_back(Lambda4); pureParDimE.push_back(1.);
    pureParVec.push_back(Cf);      pureParDimE.push_back(2.);

    


    /* additional regressors besides the main regressor */
    SymbolStringVec ssvec;
    ssvec.push_back( SymbolStringPair( &L , "L" ) );
    ssvec.push_back( SymbolStringPair( &ZP , "ZP" ) );

    return    eval_ex_lso(fpisq, /* the expression to work on */
			  x,par,aargs,deri, /* pass on the parameters from R environment */
			  mu, /* the main regressor appearing in the expression */
			  pureParVec,  /* a vector of parameters to optimize for */
			  pureParDimE,
			  ssvec   /* a vector of additional regresssor and their name in the aargs list */
			  );
  }








};
