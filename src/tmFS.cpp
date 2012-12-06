




// needed for the symbols


#include <ginac/ginac.h>

using namespace GiNaC;

#include <Rcpp.h>

#include "mpi.mq.ob.h"


#include "tmFS.h"
#include "gtilde.h"
#include "tmFS.B20.fn.h"
#include "tmFS.R.fn.h"
#include "symbols.h"
#include "eval.ex.lso.h"

namespace chifit{

  ex get_tm_FSE(ParameterMap & pm){
    ex lambda_pm=sqrt( Mpm_sq ) *  L;
    ex lambda_0=sqrt( M0_sq ) *  L;

    ex r_0 = sqrt( M0_sq / Mpm_sq );

    static ex xi_pm = Mpm_sq/ pow( 4. * Pi * f , 2 );
    static ex xi_0 = M0_sq/ pow( 4. * Pi * f , 2 );

    static ex log_Mpm_L1 = log( Mpm_sq / pow( Lambda1 , 2 ) );
    static ex log_Mpm_L2 = log( Mpm_sq / pow( Lambda2 , 2 ) );
    static ex log_Mpm_L3 = log( Mpm_sq / pow( Lambda3 , 2 ) );
    static ex log_Mpm_L4 = log( Mpm_sq / pow( Lambda4 , 2 ) );



    ex I_M_0_2_contr = xi_0/2.*gtilde1(lambda_0); // this corresponds to the ordinary XPT Xpres. but with the neutral pion mass

    ex I_M_0_4_B_0_contr = 
      - xi_0 * xi_pm 
      * gtilde1(lambda_0) 
      * (
	 -4./3.*log_Mpm_L1
	 +1./2.*log_Mpm_L3
	 +2.*log_Mpm_L4
	 + 13./18.
	 ) ;

    ex I_M_0_4_B_2_contr = 
      - 2 * pow(xi_pm,2) / lambda_pm 
      * (20./9.+8./3.*log_Mpm_L2)
      * pow(r_0,3)
      * B2k_0( lambda_0 , 1 );

    ex I_M_0_4_Rs_contr = 
      - 2 * pow(xi_pm,2) / lambda_pm
      * 2./3.
      *(
	 1. *     r_0    * tmFS_R(lambda_0,0,r_0)
	+2. * pow(r_0,2) * tmFS_R(lambda_0,1,r_0)
	-4. * pow(r_0,3) * tmFS_R(lambda_0,2,r_0)
	);



    // I_M_pm_2_contr = 0  !! -> no work needed

    ex I_M_pm_4_B_0_contr = 
      - pow( xi_pm , 2) 
      * ( - 8./3.*(  log_Mpm_L1 + log_Mpm_L2 ) + 2. * log_Mpm_L3 - 34./9.)
      * gtilde1(lambda_pm);


    ex I_M_pm_4_B_2_contr = 
      - 2 * pow(xi_pm,2) / lambda_pm 
      * (92./9.+8./3.*log_Mpm_L1 + 8. * log_Mpm_L2)
      // * r_pm^3 = 1 (r_pm = 1) 
      * B2k_0( lambda_pm , 1 );



    ex I_M_pm_4_Rs_contr = 
      - 2 * pow(xi_pm,2) / lambda_pm
      * 1./3.
      *(
	11. * tmFS_R(lambda_pm,0,1.)
	-20. * tmFS_R(lambda_pm,1,1.)
	-32. * tmFS_R(lambda_pm,2,1.)
	)
      ;


    pm.add( B );
    pm.add( f );
    pm.add( c2 );
    pm.add( Lambda1 );
    pm.add( Lambda2 );
    pm.add( Lambda3 );
    pm.add( Lambda4 );

    ex X=
        I_M_0_2_contr
      + I_M_0_4_B_0_contr
      + I_M_0_4_B_2_contr
      + I_M_0_4_Rs_contr
      + I_M_pm_4_B_0_contr
      + I_M_pm_4_B_2_contr
      + I_M_pm_4_Rs_contr
      ;





    return X;

  }


  RcppExport SEXP evalTMFSE(SEXP par,SEXP mu,SEXP L,SEXP ZP){
    Rcpp::NumericVector vpar(par);
    ParameterMap pm;

    ex X=get_tm_FSE(pm);


    const SymbolBoolVec & useMap=pm.getMap();

    exmap subsmap;

    subsmap[chifit::mu] = Rcpp::as<double>(mu);
    subsmap[chifit::L] = Rcpp::as<double>(L);
    subsmap[chifit::ZP] = Rcpp::as<double>(ZP);

    for( int p  = 0, lincount= 0; p < useMap.size() ; p++){
      if( useMap[p].second ) { 
	subsmap[useMap[p].first] = vpar[lincount++];
      }
    }



    return Rcpp::wrap( numEvalXPression( X, subsmap));


  }

};
