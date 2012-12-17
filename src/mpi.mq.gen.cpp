
#include <vector>
#include <iostream>
#include <iomanip>

using namespace std;

#include <Rcpp.h>

using namespace Rcpp;

#include "symbols.h"
#include "mpi.mq.gen.h"
#include "eval.ex.lso.h"
#include "ExpressionFactory.h"

namespace chifit {


  SEXP mpi_mq_gen(ExGenFN getEx,SEXP x, SEXP par,SEXP aargs,SEXP deri,SEXP fitZP,ExGenFNFSE fseFN) {

    ParameterMap pm;

    GiNaC::ex mpisq;

    if(fseFN == NULL) 
      mpisq = getEx(pm,true);
    else
      mpisq = getEx(pm,false)*fseFN(pm);

//    pm.print();

    /* the main parameters to optimize for */
    SymbolVec pureParVec;
    vector<double> pureParDimE;

    /* the parameters to be fitted and their energy dimension */
    /* they will be autagically selected by the expression generating function */

    const SymbolBoolVec & useMap=pm.getMap();


    for( int p  = 0 ; p < useMap.size() ; p++){
      if( useMap[p].second ) { 
	//	cout <<  useMap[p].first << endl ;
	pureParVec.push_back( useMap[p].first );
	pureParDimE.push_back( allSymbolsOrderedDimensions[p] );
      }
    }
    




    /* additional regressors besides the main regressor */
    SymbolStringVec ssvec;
    SymbolVec lsDepPar;
    ssvec.push_back( SymbolStringPair( &L , "L" ) );

    if( ! Rcpp::as<bool>(fitZP) )
      ssvec.push_back( SymbolStringPair( &ZP , "ZP" ) );
    else      
      lsDepPar.push_back(ZP);


    return    eval_ex_lso(mpisq, /* the expression to work on */
			  x,par,aargs,deri, /* pass on the parameters from R environment */
			  mu, /* the main regressor appearing in the expression */
			  pureParVec,  /* a vector of parameters to optimize for */
			  pureParDimE,
			  ssvec,   /* a vector of additional regresssor and their name in the aargs list */
			  lsDepPar);
  }

  RcppExport SEXP mpi_mq_gen_id(SEXP id,SEXP x, SEXP par,SEXP aargs,SEXP deri,SEXP fitZP) {

    ParameterMap pm;

    GiNaC::ex mpisq;

    mpisq = ExpressionFactory::getInstance()->getExpression(as<int>(id));
    pm = ExpressionFactory::getInstance()->getParameterMap(as<int>(id));


    /* the main parameters to optimize for */
    SymbolVec pureParVec;
    vector<double> pureParDimE;

    /* the parameters to be fitted and their energy dimension */
    /* they will be autagically selected by the expression generating function */

    const SymbolBoolVec & useMap=pm.getMap();


    for( int p  = 0 ; p < useMap.size() ; p++){
      if( useMap[p].second ) { 
	//	cout <<  useMap[p].first << endl ;
	pureParVec.push_back( useMap[p].first );
	pureParDimE.push_back( allSymbolsOrderedDimensions[p] );
      }
    }
    




    /* additional regressors besides the main regressor */
    SymbolStringVec ssvec;
    SymbolVec lsDepPar;
    ssvec.push_back( SymbolStringPair( &L , "L" ) );

    if( ! Rcpp::as<bool>(fitZP) )
      ssvec.push_back( SymbolStringPair( &ZP , "ZP" ) );
    else      
      lsDepPar.push_back(ZP);


    return    eval_ex_lso(mpisq, /* the expression to work on */
			  x,par,aargs,deri, /* pass on the parameters from R environment */
			  mu, /* the main regressor appearing in the expression */
			  pureParVec,  /* a vector of parameters to optimize for */
			  pureParDimE,
			  ssvec,   /* a vector of additional regresssor and their name in the aargs list */
			  lsDepPar);
  }

};
