

#ifndef TMFS_R_FN_H
#define TMFS_R_FN_H

#include <ginac/ginac.h>

using namespace GiNaC;

typedef double (*IntegrandFN)(double,void*) ;

DECLARE_FUNCTION_3P(tmFS_R);
DECLARE_FUNCTION_3P(tmFS_R_dx);
//DECLARE_FUNCTION_3P(tmFS_R_dr);

#endif
