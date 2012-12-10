

#ifndef MPI_MQ_GEN_H
#define MPI_MQ_GEN_H


#include <ginac/ginac.h>

#include "symbols.h"

namespace chifit{

  typedef GiNaC::ex(*ExGenFN)(ParameterMap &,bool);

};
#endif
