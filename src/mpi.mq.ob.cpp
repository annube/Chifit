#include <ginac/ginac.h>

#include <iostream>
#include <string>
#include <Rcpp.h>

#include <vector>

using namespace std;
using namespace GiNaC;


#include "eval.ex.h"


/* we need a namespace here because the same symbols as below are defined in mpi.mq.cpp as well */
namespace mpi_mq_ob {

  symbol B0("B_0");
  symbol f("f");
  symbol a("a");
  symbol amu_q("a\\mu_q");
  symbol Lambda3("\\lambda_3");
  symbol CMpm("C_{M^\\pm}");


  ex Mpm_sq = 2*B0*amu_q/a; 
  ex M0_sq = 2*B0*amu_q/a; 


  ex mpisq = Mpm_sq * ( 1 
			+ M0_sq / 2 / pow( 4 * Pi * f , 2 )  * log( M0_sq / pow(Lambda3,2) ) 
			+ CMpm * pow( a , 2) )  ;


  typedef map< string  , vector<double> > AargsType;


  RcppExport SEXP mpi_mq_ob(SEXP x, SEXP par,SEXP aargs) {

    /* convert SEXP's to useful Rcpp objects */
    Rcpp::NumericVector vpar(par);
    Rcpp::NumericVector vx(x);
    Rcpp::NumericVector vres(vx.size());
    Rcpp::List aargsMap ( aargs);
    Rcpp::NumericVector latSpac = aargsMap["a"] ;

    /* create a Symbol to input parameter mapping */
    exmap map;
    map[B0] = vpar[0];
    map[f] = vpar[1];
    map[Lambda3] = vpar[2];
    map[CMpm] = vpar[3];

    //void eval_ex(ex Xpress, const exmap &parValues,const XValueMap &xs,vector<double> &result);
    XValueMap xs;
    xs[ amu_q ] = &vx;
    xs[ a ] = &latSpac;

    eval_ex(mpisq,map,xs,vres);

    return vres;


  }



  template <typename T> class setSize 
  {

  private:
    int size;

  public:
    setSize(int n):size(n){}

    void operator()( T &obj ){ obj.reserve(size); }

  };


  typedef vector<symbol*>::iterator exmapIt;


  RcppExport SEXP dmpi_mq_ob(SEXP x, SEXP par,SEXP aargs) {

    Rcpp::NumericVector vpar(par);
    Rcpp::NumericVector vx(x);
    Rcpp::List aargsMap ( aargs);
    Rcpp::NumericVector latSpac = aargsMap["a"] ;


    Rcpp::NumericMatrix gradRes(vx.size(),vpar.size());


    /* create a Symbol to input parameter mapping */
    exmap map;
    map[B0] = vpar[0];
    map[f] = vpar[1];
    map[Lambda3] = vpar[2];
    map[CMpm] = vpar[3];


    //void eval_ex(ex Xpress, const exmap &parValues,const XValueMap &xs,vector<double> &result);
    XValueMap xs;
    xs[ amu_q ] = &vx;
    xs[ a ] = &latSpac;

    vector<symbol*> parvec(4);
    parvec[0] = &B0;
    parvec[1] = &f;
    parvec[2] = &Lambda3;
    parvec[3] = &CMpm;


    eval_ex_deri(mpisq,map,parvec,xs,gradRes);




    return Rcpp::wrap( gradRes );


  }

};
