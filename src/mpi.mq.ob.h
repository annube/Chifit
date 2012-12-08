


#ifndef MPI_MQ_OB_H
#define MPI_MQ_OB_H

/** 
 * function definitions for obtaining the corresponding expressions 
 */


#include <ginac/ginac.h>
#include "symbols.h"

namespace chifit {

  GiNaC::ex get_M_pm_sq_Xpression(ParameterMap &pm);
  GiNaC::ex get_M_0_sq_Xpression();
  GiNaC::ex get_f_pm_Xpression();


class Mpm_sq_ex{public: 
    Mpm_sq_ex(){} 
  operator GiNaC::ex() const {return 2*B*mu/ZP;}
};

class M0_sq_ex{public: 
    M0_sq_ex(){} 
  operator GiNaC::ex() const {return sqrt( pow(  2*B*mu/ZP + 2*c2 ,2 ) );}
};

};

#endif
