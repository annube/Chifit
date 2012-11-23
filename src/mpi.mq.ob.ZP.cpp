





#include <ginac/ginac.h>

#include <vector>

#include "eval.ex.lso.h"
#include "mpi.mq.ob.h"


using namespace GiNaC;
using namespace std;


namespace mpi_mq_ob{


RcppExport SEXP mpi_mq_ob_zp(SEXP x, SEXP par,SEXP aargs,SEXP deri) {
  static  ex mpisq = get_M_pm_sq_Xpression();
  
  /* the main parameters to optimize for */
  SymbolVec pureParVec;
  vector<double> pureParDimE;
  
  /* the parameters to be fitted and their energy dimension */
  
  pureParVec.push_back(B);       pureParDimE.push_back(1.);
  pureParVec.push_back(f);       pureParDimE.push_back(1.);
  pureParVec.push_back(c2);      pureParDimE.push_back(4.);
  pureParVec.push_back(Lambda3); pureParDimE.push_back(1.);
  pureParVec.push_back(CMpm);    pureParDimE.push_back(2.);

    


    /* additional regressors besides the main regressor */
    SymbolStringVec ssvec;
    ssvec.push_back( SymbolStringPair( &L , "L" ) );

    SymbolVec lsDepPar;
    lsDepPar.push_back(ZP);



    return    eval_ex_lso(mpisq, /* the expression to work on */
			  x,par,aargs,deri, /* pass on the parameters from R environment */
			  mu, /* the main regressor appearing in the expression */
			  pureParVec,  /* a vector of parameters to optimize for */
			  pureParDimE,
			  ssvec,   /* a vector of additional regresssor and their name in the aargs list */
			  lsDepPar
			  );
  }


};