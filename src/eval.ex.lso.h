


#ifndef EVAL_EX_LSO_H
#define EVAL_EX_LSO_H

#include <map>
#include <vector>

#include <string>
#include <utility>

#include <ginac/ginac.h>

#include <Rcpp.h>


//
typedef std::vector< GiNaC::symbol > SymbolVec;
typedef SymbolVec::iterator SymbolVecIt;

// types for additional regressors

typedef std::pair<GiNaC::symbol *,std::string> SymbolStringPair;
typedef std::vector< SymbolStringPair > SymbolStringVec;
typedef SymbolStringVec::iterator SymbolStringVecIt;


typedef std::pair<GiNaC::symbol *,Rcpp::NumericVector> SymbolVecPair;
typedef std::vector< SymbolVecPair > SymbolVecVec;
typedef SymbolVecVec::iterator SymbolVecVecIt;

typedef std::vector<GiNaC::ex> ExVec;
typedef ExVec::iterator ExVecIt;


namespace chifit {

double numEvalXPression(GiNaC::ex Xpression,GiNaC::exmap substitutions);


SEXP eval_ex_lso(
		 const GiNaC::ex &Xpression,SEXP x,SEXP par, SEXP aargs,SEXP deri,
		 GiNaC::symbol mainRegressor,
		 SymbolVec svec,
		 std::vector<double> pureParDimE,
		 SymbolStringVec addRegr,
		 SymbolVec lsDepPar = SymbolVec(0)
		 ) ;
};



#endif
