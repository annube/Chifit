#include <ginac/ginac.h>

#include <iostream>
#include <string>
#include <Rcpp.h>

#include <vector>

using namespace std;

using namespace GiNaC;


#include "gtilde.h"
#include "eval.ex.h"

namespace mpi_mq_lso{

  symbol aB0("aB_0");
  symbol af("af");
  symbol amu_q("a\\mu_q");
  symbol aLambda3("a\\Lambda_3");
  symbol aLambda4("a\\Lambda_4");
  symbol L("L");
  symbol ZP("Z_P");

  ex chi_mu_asq = 2.*aB0*amu_q/ZP; 
  ex xi_ll = chi_mu_asq / pow( 4. * Pi * af , 2. );



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
    static ex X =  chi_mu_asq * ( 1.
					   + xi_ll  * log( chi_mu_asq / pow( aLambda3,2.) ) 
					   ) * pow(1. + 0.5 * xi_ll * gtilde1( sqrt(
										 chi_mu_asq ) *  L
									       ) ,2) ;
    return X;
  }

  RcppExport SEXP mpi_mq_lso(SEXP x, SEXP par,SEXP aargs,SEXP deri) {


  // pion mass
    static  ex ampisq = getampisqXpression();

    /* convert SEXP's to useful Rcpp objects */
    Rcpp::NumericVector vpar(par);
    Rcpp::NumericVector vx(x);

    Rcpp::NumericVector vres(vx.size());
    Rcpp::NumericMatrix gradRes(vx.size(),vpar.size());
    Rcpp::NumericMatrix gradRaw(vx.size(),3);

    bool calcDeri = Rcpp::as<bool>(deri);

    Rcpp::List aargsMap ( aargs);
    Rcpp::NumericVector numLs = aargsMap["numLs"] ;
    Rcpp::NumericVector latSpacIndex = aargsMap["lsIndex"] ;
    Rcpp::NumericVector Loa = aargsMap["L"] ;
    Rcpp::NumericVector zpVal = aargsMap["ZP"] ;
    

    Rcpp::NumericVector af0Vals(vx.size());
    Rcpp::NumericVector aB0Vals(vx.size());
    Rcpp::NumericVector aLambda3Vals(vx.size());

    vector<double> lsRatios(numLs[0]);

    for(int j = 0 ; j < numLs[0] ; j++)
      lsRatios[j] = vpar[ 1 + j ]/vpar[1];


    for(int i = 0 ; i<vx.size() ; i++){
      af0Vals[i] = vpar[latSpacIndex[i] + 0 /* R's indices start from 1 */ ];
      aB0Vals[i] = vpar[0] * lsRatios[ latSpacIndex[i] - 1 ];
      aLambda3Vals[i] = vpar[1 + numLs[0]] * lsRatios[ latSpacIndex[i] - 1 ];
    }

    /* create a Symbol to input parameter mapping */
    exmap map;



    //void eval_ex(ex Xpress, const exmap &parValues,const XValueMap &xs,vector<double> &result);
    XValueMap xs;
    xs[ af ] = &af0Vals;
    xs[ aB0 ] = &aB0Vals;
    xs[ aLambda3 ] = &aLambda3Vals;
    xs[ amu_q ] = &vx;
    xs[L] = &Loa;
    xs[ZP] = &zpVal;

     vector<symbol*> smap(3);
     smap[0] = &aB0;
     smap[1] = &af;
     smap[2] = &aLambda3;  // <-- this is correct



    if(calcDeri){
      eval_ex_deri(ampisq,map,smap,xs,gradRaw);

      for( int ix = 0 ; ix < vx.size() ; ix++){
	for( int apar = 0 ; apar < 3 ; apar++){
	  cout << ix << " " << apar << endl;
	  switch ( apar ){
	  case 0: gradRes(ix,0) = gradRaw(ix,0); break;
	  case 1: gradRes(ix,latSpacIndex[ix]) = gradRaw(ix,1); break;
	  case 2: gradRes(ix,1+numLs[0]) = gradRaw(ix,2); break;

	  }
	}

      }

      return gradRes;
    } else {
      eval_ex(ampisq,map,xs,vres);
      return vres;
    }



  }



//   RcppExport SEXP dmpi_mq_lso(SEXP x, SEXP par,SEXP aargs) {

//   // pion mass
//     static  ex ampisq = getampisqXpression();

//     Rcpp::NumericVector vpar(par);
//     Rcpp::NumericVector vx(x);
//     Rcpp::List aargsMap ( aargs);
//     Rcpp::NumericVector latSpac = aargsMap["a"] ;
//     Rcpp::NumericVector Loa = aargsMap["L"] ;
//     Rcpp::NumericVector zpVal = aargsMap["ZP"] ;


//     Rcpp::NumericMatrix gradRes(vx.size(),vpar.size());


//     /* create a Symbol to input parameter mapping */
//     exmap map;
//     map[B0] = vpar[0];
//     map[f] = vpar[1];
//     map[Lambda3] = vpar[2];

//     vector<symbol*> smap(3);
//     smap[0] = &B0;
//     smap[1] = &f;
//     smap[2] = &Lambda3;  // <-- this is correct

//     //void eval_ex(ex Xpress, const exmap &parValues,const XValueMap &xs,vector<double> &result);
//     XValueMap xs;
//     xs[ amu_q ] = &vx;
//     xs[ a ] = &latSpac;
//     xs[ L ] = &Loa;
//     xs[ ZP ] = &zpVal;


//     eval_ex_deri(ampisq,map,smap,xs,gradRes);





//     return Rcpp::wrap( gradRes );


//   }


//   /**************************
//    *
//    *
//    *  PION DECAY CONSTANT
//    *
//    *
//    *
//    **************************/




//   static ex getafpiXpression(){
//     static ex X = a * f * ( 1 
// 			    - 2 *  xi_ll  * log( chi_mu / pow( Lambda4,2) ) 
// 			    )
//       *
//       ( 1 - 2  * xi_ll * gtilde1( sqrt(  chi_mu ) *  L * a   ) ) ;
//     return X;
//   }


//   RcppExport SEXP fpi_mq_lso(SEXP x, SEXP par,SEXP aargs) {

//     /* convert SEXP's to useful Rcpp objects */
//     Rcpp::NumericVector vpar(par);
//     Rcpp::NumericVector vx(x);
//     Rcpp::NumericVector vres(vx.size());
//     Rcpp::List aargsMap ( aargs);
//     Rcpp::NumericVector latSpac = aargsMap["a"] ;
//     Rcpp::NumericVector Loa = aargsMap["L"] ;
//     Rcpp::NumericVector zpVal = aargsMap["ZP"] ;

//     // pion decay
//     static ex afpi = getafpiXpression();


//     /* create a Symbol to input parameter mapping */
//     exmap map;
//     map[B0] = vpar[0];
//     map[f] = vpar[1];
//     map[Lambda4] = vpar[2];


//     //void eval_ex(ex Xpress, const exmap &parValues,const XValueMap &xs,vector<double> &result);
//     XValueMap xs;
//     xs[ amu_q ] = &vx;
//     xs[ a ] = &latSpac;
//     xs[ L ] = &Loa;
//     xs[ ZP ] = &zpVal;


//     eval_ex(afpi,map,xs,vres);


//     return vres;


//   }



//   RcppExport SEXP dfpi_mq_lso(SEXP x, SEXP par,SEXP aargs) {

//     Rcpp::NumericVector vpar(par);
//     Rcpp::NumericVector vx(x);
//     Rcpp::List aargsMap ( aargs);
//     Rcpp::NumericVector latSpac = aargsMap["a"] ;
//     Rcpp::NumericVector Loa = aargsMap["L"] ;
//     Rcpp::NumericVector zpVal = aargsMap["ZP"] ;


//     Rcpp::NumericMatrix gradRes(vx.size(),vpar.size());


//     // pion decay
//     static ex afpi = getafpiXpression();

//     /* create a Symbol to input parameter mapping */
//     exmap map;
//     map[B0] = vpar[0];
//     map[f] = vpar[1];
//     map[Lambda4] = vpar[2];

//     vector<symbol*> smap(3);
//     smap[0] = &B0;
//     smap[1] = &f;
//     smap[2] = &Lambda4;

//     //void eval_ex(ex Xpress, const exmap &parValues,const XValueMap &xs,vector<double> &result);
//     XValueMap xs;
//     xs[ amu_q ] = &vx;
//     xs[ a ] = &latSpac;
//     xs[ L ] = &Loa;
//     xs[ ZP ] = &zpVal;


//     eval_ex_deri(afpi,map,smap,xs,gradRes);





//     return Rcpp::wrap( gradRes );


//   }

};
