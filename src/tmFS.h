

#ifndef TMFS_H
#define TMFS_H

#include <ginac/ginac.h>

#include "symbols.h"

namespace chifit{

  GiNaC::ex get_tm_FSE_Rpm(ParameterMap & pm,bool simplified = true);
  GiNaC::ex get_tm_FSE_Mpm_sq ( ParameterMap & pm);
  GiNaC::ex get_tm_FSE_RM0(ParameterMap & pm);
  GiNaC::ex get_tm_FSE_Rfpm(ParameterMap & pm, bool simplified = true);

}


#endif
