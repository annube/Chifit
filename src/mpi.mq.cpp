#include <ginac/ginac.h>

#include <iostream>
#include <string>
#include <Rcpp.h>

#include <vector>

using namespace std;

using namespace GiNaC;


#include "eval.ex.h"

namespace mpi_mq{

  symbol B0("B_0");
  symbol f("f");
  symbol a("a");
  symbol amu_q("a\\mu_q");
  symbol Lambda3("\\lambda_3");

  ex chi_mu = 2*B0*amu_q/a; 


  ex ampisq = pow(a , 2) * chi_mu * ( 1 
				      + pow(a,  2) * chi_mu / pow( 4 * Pi * a * f , 2 )  * log( pow(a,2)*chi_mu / pow( a * Lambda3,2) ) 
				      )  ;


  typedef map< string  , vector<double> > AargsType;


  RcppExport SEXP mpi_mq(SEXP x, SEXP par,SEXP aargs) {

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


    //void eval_ex(ex Xpress, const exmap &parValues,const XValueMap &xs,vector<double> &result);
    XValueMap xs;
    xs[ amu_q ] = &vx;
    xs[ a ] = &latSpac;


    eval_ex(ampisq,map,xs,vres);


    return vres;


  }



  typedef vector<symbol*>::iterator exmapIt;


  RcppExport SEXP dmpi_mq(SEXP x, SEXP par,SEXP aargs) {

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

    vector<symbol*> smap(3);
    smap[0] = &B0;
    smap[1] = &f;
    smap[2] = &Lambda3;

    //void eval_ex(ex Xpress, const exmap &parValues,const XValueMap &xs,vector<double> &result);
    XValueMap xs;
    xs[ amu_q ] = &vx;
    xs[ a ] = &latSpac;


    eval_ex_deri(ampisq,map,smap,xs,gradRes);





    return Rcpp::wrap( gradRes );


  }

};
