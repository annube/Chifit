




// needed for the symbols


#include <ginac/ginac.h>

using namespace GiNaC;


#include "mpi.mq.ob.h"

using namespace mpi_mq_ob;

#include "tmFS.h"

#include "gtilde.h"


namespace tmFS{

  symbol Lambda1("\\Lambda_1");
  symbol Lambda2("\\Lambda_2");



  ex get_tm_FSE(){
    ex lambda_pm=sqrt( Mpm_sq ) *  L;
    ex lambda_0=sqrt( M0_sq ) *  L;

    static ex xi_pm = Mpm_sq/ pow( 4. * Pi * f , 2 );
    static ex xi_0 = M0_sq/ pow( 4. * Pi * f , 2 );

    static ex log_Mpm_L1 = log( Mpm_sq / pow( Lambda1 , 2 ) );
    static ex log_Mpm_L2 = log( Mpm_sq / pow( Lambda2 , 2 ) );
    static ex log_Mpm_L3 = log( Mpm_sq / pow( Lambda3 , 2 ) );
    static ex log_Mpm_L4 = log( Mpm_sq / pow( Lambda4 , 2 ) );

    ex I_M_0_2_contr = xi_0/2.*gtilde1(lambda_0);
    ex I_M_0_4_B_0_contr = 
      - xi_0 * xi_pm 
      * gtilde1(lambda_0) 
      * (
	 -4./3.*log_Mpm_L1
	 +1./2.*log_Mpm_L3
	 +2.*log_Mpm_L4
	 + 13./18.
	 )
      * gtilde1(lambda_0)
      ;


    return 0;

  }

};
