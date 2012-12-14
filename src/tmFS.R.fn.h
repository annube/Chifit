

#ifndef TMFS_R_FN_H
#define TMFS_R_FN_H

#include <ginac/ginac.h>

using namespace GiNaC;

typedef double (*IntegrandFN)(double,void*) ;

typedef struct R_Integrand_Params_ {
  int k,n;
  double x,r;
} R_Integrand_Params;


double R_integrand(double y, void *params);


double eval_tmFS_R(double x,int k,double r,IntegrandFN integrand=&R_integrand);

DECLARE_FUNCTION_3P(tmFS_R);
DECLARE_FUNCTION_3P(tmFS_R_dx);
DECLARE_FUNCTION_3P(tmFS_R_dr);

#endif
