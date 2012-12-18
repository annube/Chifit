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
#include "mpi.mq.gen.h"
#include "tmFS.h"





/* we need a namespace here because the same symbols as below are defined in mpi.mq.cpp as well */
namespace chifit {

  Mpm_sq_ex Mpm_sq;
  M0_sq_ex M0_sq; 


  ex get_M_pm_sq_ob_Xpression(ParameterMap & pm , bool withFS /* =true */ ){
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


  ex get_M_0_sq_ob_Xpression(ParameterMap &pm,bool withFS){

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

    pm.add( B );
    pm.add( f );
    pm.add( c2 );
    pm.add( Lambda3 );
    pm.add( Xi3 );
    pm.add( CM0 );

    if( withFS )
      return X_FSE;
    else
      return X;

  }


  ex get_f_pm_ob_Xpression(ParameterMap &pm,bool withFS){
    static ex xi_pm = Mpm_sq/ pow( 4. * Pi * f , 2 );
    static ex xi_0  = M0_sq/ pow( 4. * Pi * f , 2 );
    static ex fse_log_corr_Mpm = gtilde1( sqrt( Mpm_sq ) *  L );
    static ex fse_log_corr_M0  = gtilde1( sqrt( M0_sq ) *  L );
    static ex log_M0_L4 = log( M0_sq / pow( Lambda4 , 2 )  );
    static ex log_Mpm_L4 = log( Mpm_sq / pow( Lambda4 , 2 ) );

    static ex X_FSE =  f * ( 1. - ( xi_pm * ( log_Mpm_L4 + fse_log_corr_Mpm ) + xi_0 * (  log_M0_L4  + fse_log_corr_M0 ) ) + Cf ) ;

    static ex X =  f * ( 1. - ( xi_pm * log_Mpm_L4 + xi_0 * log_M0_L4 ) + Cf ) ;

    pm.add( B );
    pm.add( f );
    pm.add( c2 );
    pm.add( Lambda4 );
    pm.add( Cf );

    if( withFS )
      return X_FSE;
    else
      return X;

  }







};
