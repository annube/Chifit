

#ifndef TMFS_H
#define TMFS_H

#include <ginac/ginac.h>

#include "symbols.h"

namespace chifit{

  GiNaC::ex get_tm_FSE(ParameterMap & pm);
  GiNaC::ex get_tm_FSE_Mpm_sq ( ParameterMap & pm);

}


#endif
