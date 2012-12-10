

#ifndef MPI_MQ_GEN_H
#define MPI_MQ_GEN_H


#include <ginac/ginac.h>

#include "symbols.h"

namespace chifit{

  typedef GiNaC::ex(*ExGenFN)(ParameterMap &,bool);
  typedef GiNaC::ex(*ExGenFNFSE)(ParameterMap &);

  SEXP mpi_mq_gen(ExGenFN getEx,SEXP x, SEXP par,SEXP aargs,SEXP deri,SEXP fitZP,ExGenFNFSE fseFN=NULL);

};
#endif
