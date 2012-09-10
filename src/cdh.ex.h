

#ifndef CDH_EX_H
#define CDH_EX_H


#include <ginac/ginac.h>

using namespace GiNaC;

ex cdh_mpi_ex(
	   const ex &mpi,
	   const ex &f0,
	   const ex &L,
	   const ex &Lambda1,
	   const ex &Lambda2,
	   const ex &Lambda3,
	   const ex &Lambda4,
	   int maxorder = 20
	   ) ;



ex cdh_fpi_ex(
	   const ex &mpisq,
	   const ex &f0,
	   const ex &L,
	   const ex &Lambda1,
	   const ex &Lambda2,
	   const ex &Lambda3,
	   const ex &Lambda4,
	   int maxorder = 20
	   );



#endif
