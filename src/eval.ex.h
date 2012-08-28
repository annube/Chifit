
#include <Rcpp.h>


typedef map< ex , Rcpp::NumericVector*,ex_is_less > XValueMap;

typedef map< ex , Rcpp::NumericVector*, ex_is_less >::const_iterator XValueMapIt;


typedef vector<symbol*> parList;
typedef vector<symbol*>::iterator parListIt;

void eval_ex(ex Xpress, const exmap &parValues,const XValueMap &xs,Rcpp::NumericVector &result);
void eval_ex_deri(ex Xpress, const exmap &parValues,parList &smap,const XValueMap &xs,Rcpp::NumericMatrix &result) ;

