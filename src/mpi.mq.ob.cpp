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
#include "mpi.mq.ob.h"
#include "tmFS.h"





/* we need a namespace here because the same symbols as below are defined in mpi.mq.cpp as well */
namespace chifit {

// class Mpm_sq_ex{public: 
//     Mpm_sq_ex(){} 
//   operator GiNaC::ex() const {return 2*B*mu/ZP;}
// };

// class M0_sq_ex{public: 
//     M0_sq_ex(){} 
//   operator GiNaC::ex() const {return sqrt( pow(  2*B*mu/ZP + 2*c2 ,2 ) );}
// };



  Mpm_sq_ex Mpm_sq;
  M0_sq_ex M0_sq; 


  ex get_M_pm_sq_Xpression(ParameterMap & pm , bool withFS /* =true */ ){
    static ex xi_0  = M0_sq/ pow( 4. * Pi * f , 2 );
    static ex fse_log_corr_M0  = gtilde1( sqrt( M0_sq ) *  L );
    static ex log_M0_L3 = log( M0_sq / pow( Lambda3 , 2 )  )  ;

    static ex X_FSE = (  Mpm_sq * ( 1. + xi_0 * ( log_M0_L3 + fse_log_corr_M0 )+ CMpm ) ) ; 

    static ex X = (  Mpm_sq * ( 1. + xi_0 * log_M0_L3 + CMpm ) ) ; 

    pm.add( B );
    pm.add( f );
    pm.add( c2 );
    pm.add( Lambda3 );
    pm.add( CMpm );

    if(withFS)
      return X_FSE;
    else
      return X;


  }


  ex get_M_0_sq_Xpression(ParameterMap &pm,bool withFS){

    static ex xi_0  = M0_sq/ pow( 4. * Pi * f , 2 );
    static ex fse_log_corr_M0  = gtilde1( sqrt( M0_sq ) *  L );
    static ex log_M0_L3 = log( M0_sq / pow( Lambda3 , 2 )  ) ;

    static ex xi_pm  = Mpm_sq/ pow( 4. * Pi * f , 2 );
    static ex fse_log_corr_Mpm  = gtilde1( sqrt( Mpm_sq ) *  L );
    static ex log_Mpm_L3 = log( Mpm_sq / pow( Lambda3 , 2 )  )  ;

    static ex log_M0_X3 = log( M0_sq / pow( Xi3 , 2 )  )  ;

    static ex X_FSE =   Mpm_sq * ( 1. +2.*xi_pm * ( log_Mpm_L3 + fse_log_corr_Mpm ) - xi_0 * ( log_M0_L3 +  + fse_log_corr_M0 )) 
      + 2 * c2 *( 1 - 4. * xi_0 * ( log_M0_X3 + fse_log_corr_M0) + CM0); 

    static ex X =   Mpm_sq * ( 1. +2.*xi_pm * log_Mpm_L3 - xi_0 * log_M0_L3 ) 
      + 2 * c2 *( 1 - 4. * xi_0 * log_M0_X3 + CM0); 

    if( withFS )
      return X_FSE;
    else
      return X;

  }


  ex get_f_pm_Xpression(ParameterMap &pm,bool withFS){
    static ex xi_pm = Mpm_sq/ pow( 4. * Pi * f , 2 );
    static ex xi_0  = M0_sq/ pow( 4. * Pi * f , 2 );
    static ex fse_log_corr_Mpm = gtilde1( sqrt( Mpm_sq ) *  L );
    static ex fse_log_corr_M0  = gtilde1( sqrt( M0_sq ) *  L );
    static ex log_M0_L4 = log( M0_sq / pow( Lambda4 , 2 )  );
    static ex log_Mpm_L4 = log( Mpm_sq / pow( Lambda4 , 2 ) );

    static ex X_FSE =  f * ( 1. - ( xi_pm * ( log_Mpm_L4 + fse_log_corr_Mpm ) + xi_0 * (  log_M0_L4  + fse_log_corr_M0 ) ) + Cf ) ;

    static ex X =  f * ( 1. - ( xi_pm * log_Mpm_L4 + xi_0 * log_M0_L4 ) + Cf ) ;

    if( withFS )
      return X_FSE;
    else
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

  RcppExport SEXP mpi_mq_ob(SEXP x, SEXP par,SEXP aargs,SEXP deri,SEXP FSE,SEXP fitZP) {

    ParameterMap pm;

    ex mpisq = get_M_pm_sq_Xpression(pm);

    //    pm.print();

    /* the main parameters to optimize for */
    SymbolVec pureParVec;
    vector<double> pureParDimE;

    /* the parameters to be fitted and their energy dimension */
    /* they will be autagically selected by the expression generating function */

    const SymbolBoolVec & useMap=pm.getMap();


    for( int p  = 0 ; p < useMap.size() ; p++){
      if( useMap[p].second ) { 
	pureParVec.push_back( useMap[p].first );
	pureParDimE.push_back( allSymbolsOrderedDimensions[p] );
      }
    }
    




    /* additional regressors besides the main regressor */
    SymbolStringVec ssvec;
    SymbolVec lsDepPar;
    ssvec.push_back( SymbolStringPair( &L , "L" ) );

    if( ! Rcpp::as<bool>(fitZP) )
      ssvec.push_back( SymbolStringPair( &ZP , "ZP" ) );
    else      
      lsDepPar.push_back(ZP);


    return    eval_ex_lso(mpisq, /* the expression to work on */
			  x,par,aargs,deri, /* pass on the parameters from R environment */
			  mu, /* the main regressor appearing in the expression */
			  pureParVec,  /* a vector of parameters to optimize for */
			  pureParDimE,
			  ssvec,   /* a vector of additional regresssor and their name in the aargs list */
			  lsDepPar);
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
    ParameterMap pm;
    static  ex mpisq = get_M_0_sq_Xpression(pm);

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
    ParameterMap pm;
    static  ex fpisq = get_f_pm_Xpression(pm);

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
