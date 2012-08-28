#include <ginac/ginac.h>

#include <iostream>
#include <string>
#include <Rcpp.h>

#include <vector>

using namespace std;

using namespace GiNaC;


#include "gtilde.h"
#include "eval.ex.h"

namespace mpi_mq{

  symbol B0("B_0");
  symbol f("f");
  symbol a("a");
  symbol amu_q("a\\mu_q");
  symbol Lambda3("\\Lambda_3");
  symbol Lambda4("\\Lambda_4");
  symbol L("L");
  symbol ZP("Z_P");

  ex chi_mu = 2.*B0*amu_q/ZP/a; 
  ex xi_ll = chi_mu / pow( 4. * Pi * f , 2. );



  typedef map< string  , vector<double> > AargsType;




  /**************************
   *
   *
   *        PION MASS
   *
   *
   *
   **************************/




  static ex getampisqXpression(){
    static ex X = pow(a , 2.) * chi_mu * ( 1.
					   + xi_ll  * log( chi_mu / pow( Lambda3,2.) ) 
					   ) * pow(1. + 0.5 * xi_ll * gtilde1( sqrt(
										 chi_mu ) *  L * a
									       ) ,2) ;
    return X;
  }

  RcppExport SEXP mpi_mq(SEXP x, SEXP par,SEXP aargs) {


  // pion mass
    static  ex ampisq = getampisqXpression();

    /* convert SEXP's to useful Rcpp objects */
    Rcpp::NumericVector vpar(par);
    Rcpp::NumericVector vx(x);
    Rcpp::NumericVector vres(vx.size());
    Rcpp::List aargsMap ( aargs);
    Rcpp::NumericVector latSpac = aargsMap["a"] ;
    Rcpp::NumericVector Loa = aargsMap["L"] ;
    Rcpp::NumericVector zpVal = aargsMap["ZP"] ;
    

    /* create a Symbol to input parameter mapping */
    exmap map;
    map[B0] = vpar[0];
    map[f] = vpar[1];
    map[Lambda3] = vpar[2];



    //void eval_ex(ex Xpress, const exmap &parValues,const XValueMap &xs,vector<double> &result);
    XValueMap xs;
    xs[ amu_q ] = &vx;
    xs[ a ] = &latSpac;
    xs[L] = &Loa;
    xs[ZP] = &zpVal;


    eval_ex(ampisq,map,xs,vres);


    return vres;


  }



  RcppExport SEXP dmpi_mq(SEXP x, SEXP par,SEXP aargs) {

  // pion mass
    static  ex ampisq = getampisqXpression();

    Rcpp::NumericVector vpar(par);
    Rcpp::NumericVector vx(x);
    Rcpp::List aargsMap ( aargs);
    Rcpp::NumericVector latSpac = aargsMap["a"] ;
    Rcpp::NumericVector Loa = aargsMap["L"] ;
    Rcpp::NumericVector zpVal = aargsMap["ZP"] ;


    Rcpp::NumericMatrix gradRes(vx.size(),vpar.size());


    /* create a Symbol to input parameter mapping */
    exmap map;
    map[B0] = vpar[0];
    map[f] = vpar[1];
    map[Lambda3] = vpar[2];

    vector<symbol*> smap(3);
    smap[0] = &B0;
    smap[1] = &f;
    smap[2] = &Lambda3;  // <-- this is correct

    //void eval_ex(ex Xpress, const exmap &parValues,const XValueMap &xs,vector<double> &result);
    XValueMap xs;
    xs[ amu_q ] = &vx;
    xs[ a ] = &latSpac;
    xs[ L ] = &Loa;
    xs[ ZP ] = &zpVal;


    eval_ex_deri(ampisq,map,smap,xs,gradRes);





    return Rcpp::wrap( gradRes );


  }


  /**************************
   *
   *
   *  PION DECAY CONSTANT
   *
   *
   *
   **************************/




  static ex getafpiXpression(){
    static ex X = a * f * ( 1 
			    - 2 *  xi_ll  * log( chi_mu / pow( Lambda4,2) ) 
			    )
      *
      ( 1 - 2  * xi_ll * gtilde1( sqrt(  chi_mu ) *  L * a   ) ) ;
    return X;
  }


  RcppExport SEXP fpi_mq(SEXP x, SEXP par,SEXP aargs) {

    /* convert SEXP's to useful Rcpp objects */
    Rcpp::NumericVector vpar(par);
    Rcpp::NumericVector vx(x);
    Rcpp::NumericVector vres(vx.size());
    Rcpp::List aargsMap ( aargs);
    Rcpp::NumericVector latSpac = aargsMap["a"] ;
    Rcpp::NumericVector Loa = aargsMap["L"] ;
    Rcpp::NumericVector zpVal = aargsMap["ZP"] ;

    // pion decay
    static ex afpi = getafpiXpression();


    /* create a Symbol to input parameter mapping */
    exmap map;
    map[B0] = vpar[0];
    map[f] = vpar[1];
    map[Lambda4] = vpar[2];


    //void eval_ex(ex Xpress, const exmap &parValues,const XValueMap &xs,vector<double> &result);
    XValueMap xs;
    xs[ amu_q ] = &vx;
    xs[ a ] = &latSpac;
    xs[ L ] = &Loa;
    xs[ ZP ] = &zpVal;


    eval_ex(afpi,map,xs,vres);


    return vres;


  }



  RcppExport SEXP dfpi_mq(SEXP x, SEXP par,SEXP aargs) {

    Rcpp::NumericVector vpar(par);
    Rcpp::NumericVector vx(x);
    Rcpp::List aargsMap ( aargs);
    Rcpp::NumericVector latSpac = aargsMap["a"] ;
    Rcpp::NumericVector Loa = aargsMap["L"] ;
    Rcpp::NumericVector zpVal = aargsMap["ZP"] ;


    Rcpp::NumericMatrix gradRes(vx.size(),vpar.size());


    // pion decay
    static ex afpi = getafpiXpression();

    /* create a Symbol to input parameter mapping */
    exmap map;
    map[B0] = vpar[0];
    map[f] = vpar[1];
    map[Lambda4] = vpar[2];

    vector<symbol*> smap(3);
    smap[0] = &B0;
    smap[1] = &f;
    smap[2] = &Lambda4;

    //void eval_ex(ex Xpress, const exmap &parValues,const XValueMap &xs,vector<double> &result);
    XValueMap xs;
    xs[ amu_q ] = &vx;
    xs[ a ] = &latSpac;
    xs[ L ] = &Loa;
    xs[ ZP ] = &zpVal;


    eval_ex_deri(afpi,map,smap,xs,gradRes);





    return Rcpp::wrap( gradRes );


  }

};
