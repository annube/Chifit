


#include <math.h>
#include <ginac/ginac.h>

using namespace GiNaC;


#include <vector>

#include "cdh.ex.h"
#include "bessel.h"
#include "gtilde.partitions.h"

ex cdh_mpi_ex(
	   const ex &mpisq,
	   const ex &f0,
	   const ex &L,
	   const ex &Lambda1,
	   const ex &Lambda2,
	   const ex &Lambda3,
	   const ex &Lambda4,
	   int maxorder
	   ) {


  std::vector<ex> exps;

  exps.push_back(1.);

  for(int i = 1 ; i <= maxorder ; i++){

    ex ln = ::sqrt( (double) i ) *  sqrt( mpisq ) * L;
    ex xi= mpisq/ pow( 4 * Pi * f0 , 2 );
    ex I_m_2 = -2. * bessel_k_fn( ln ,1 );

    ex l1 = log( mpisq / pow(Lambda1,2) );
    ex l2 = log( mpisq / pow(Lambda2,2) );
    ex l3 = log( mpisq / pow(Lambda3,2) );
    ex l4 = log( mpisq / pow(Lambda4,2) );

    ex I_m_4 = (
		101./9. - 13./3. * Pi 
		+ 8. * l1 
		+ 16./3. * l2 
		- 5. * l3 
		- 4 *l4 
		) * bessel_k_fn(ln,1) 
      + (
	 - 238. / 9. + 61. / 6. * Pi
	 -  16. / 3. * l1 
	 -  64. / 3. * l2 
	 ) * bessel_k_fn(ln,2)/ln
      ;

    exps.push_back( - (double) getPartition( i ) / 2. / ln * xi * ( I_m_2 + xi * I_m_4 ) );

  }

  return add(exps);
}



ex cdh_fpi_ex(
	   const ex &mpisq,
	   const ex &f0,
	   const ex &L,
	   const ex &Lambda1,
	   const ex &Lambda2,
	   const ex &Lambda3,
	   const ex &Lambda4,
	   int maxorder
	   ) {


  std::vector<ex> exps;

  exps.push_back(1.);

  for(int i = 1 ; i <= maxorder ; i++){

    ex ln = ::sqrt( (double) i ) *  sqrt( mpisq ) * L;
    ex xi= mpisq/ pow( 4 * Pi * f0 , 2 );
    ex I_m_2 = -4. * bessel_k_fn( ln ,1 );

    ex l1 = log( mpisq / pow(Lambda1,2) );
    ex l2 = log( mpisq / pow(Lambda2,2) );
    ex l3 = log( mpisq / pow(Lambda3,2) );
    ex l4 = log( mpisq / pow(Lambda4,2) );

    ex I_m_4 = (
		29./18. - 29./12. * Pi 
		+ 4. * l1 
		+ 8./3. * l2 
		- 6. *l4 
		) * bessel_k_fn(ln,1) 
      + (
	 - 307. / 9. + 391. / 24. * Pi
	 -  16. / 3. * l1 
	 -  64. / 3. * l2 
	 ) * bessel_k_fn(ln,2)/ln
      ;

    exps.push_back( + (double) getPartition( i ) / ln * xi * ( I_m_2 + xi * I_m_4 ) );

  }

  return add(exps);
}




#include <Rcpp.h>


RcppExport SEXP eval_cdh_mpi_ex( SEXP x){

  Rcpp::NumericVector par(x);

  symbol mpisq,f0,L,L1,L2,L3,L4;


  exmap subsmap;
  subsmap[mpisq] = par[0];
  subsmap[f0] = par[1];
  subsmap[L] = par[2];
  subsmap[L1] = par[3];
  subsmap[L2] = par[4];
  subsmap[L3] = par[5];
  subsmap[L4] = par[6];

  return Rcpp::wrap(ex_to<numeric>( cdh_mpi_ex(
					    mpisq,
					    f0,
					    L,
					    L1,
					    L2,
					    L3,
					    L4
					    ).subs(
							    subsmap
							    ).evalf()).to_double() ) ;
}

