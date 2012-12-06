



/** 
 * function definitions for obtaining the corresponding expressions 
 */


#include <ginac/ginac.h>
#include "symbols.h"

namespace chifit {

  extern GiNaC::ex Mpm_sq;
  extern GiNaC::ex M0_sq;

  GiNaC::ex get_M_pm_sq_Xpression(ParameterMap &pm);
  GiNaC::ex get_M_0_sq_Xpression();
  GiNaC::ex get_f_pm_Xpression();

};
