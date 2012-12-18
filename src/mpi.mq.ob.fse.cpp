

#include <ginac/ginac.h>


using namespace GiNaC;


#include "symbols.h"
#include "mpi.mq.ob.h"
#include "gtilde.h"



/* we need a namespace here because the same symbols as below are defined in mpi.mq.cpp as well */
namespace chifit {



  ex get_M_pm_sq_ob_FSE_Xpression(ParameterMap & pm  ){

    static ex X_FSE =   Mpm_sq *  M0_sq / pow( 4. * Pi * f , 2 ) * gtilde1( sqrt( M0_sq ) *  L ) ; 


    pm.add( B );
    pm.add( f );
    pm.add( c2 );

    return X_FSE;


  }

  ex get_M_0_sq_ob_FSE_Xpression(ParameterMap & pm  ){
     static ex xi_0  = M0_sq/ pow( 4. * Pi * f , 2 );
     static ex fse_log_corr_M0  = gtilde1( sqrt( M0_sq ) *  L );

     static ex xi_pm  = Mpm_sq/ pow( 4. * Pi * f , 2 );
     static ex fse_log_corr_Mpm  = gtilde1( sqrt( Mpm_sq ) *  L );

     static ex X_FSE =   Mpm_sq * ( 2.*xi_pm *  fse_log_corr_Mpm  - xi_0 * fse_log_corr_M0 ) 
       + 2 * c2 *( - 4. * xi_0 *  fse_log_corr_M0 ); 

    pm.add( B );
    pm.add( f );
    pm.add( c2 );

     return X_FSE;
  }

  ex get_f_pm_ob_FSE_Xpression(ParameterMap & pm  ){
    static ex xi_pm = Mpm_sq/ pow( 4. * Pi * f , 2 );
    static ex xi_0  = M0_sq/ pow( 4. * Pi * f , 2 );
    static ex fse_log_corr_Mpm = gtilde1( sqrt( Mpm_sq ) *  L );
    static ex fse_log_corr_M0  = gtilde1( sqrt( M0_sq ) *  L );

    static ex X_FSE =  f * (  - ( xi_pm * fse_log_corr_Mpm  + xi_0 *  fse_log_corr_M0  ) ) ;

    pm.add( B );
    pm.add( f );
    pm.add( c2 );

    return X_FSE;
  }





};
