
#include <iostream>

using namespace std;


#include <ginac/ginac.h>

using namespace GiNaC;

#include <Rcpp.h>

using namespace Rcpp;


#include "ExpressionFactory.h"

namespace chifit {


  
  ExpressionFactory *ExpressionFactory::instance = NULL;


  ExpressionFactory::ExpressionFactory():nextId(0){}


  ExpressionFactory *ExpressionFactory::getInstance(){
    if( ExpressionFactory::instance == NULL ){
      cout << "Generating new ExpressionFactory instance .." << endl;
      ExpressionFactory::instance = new ExpressionFactory;
    }
    return ExpressionFactory::instance;
  }

  int ExpressionFactory::generateXPression(const string &quant,const string &FSE){
    return -1;
  }




  RcppExport SEXP GenerateChifitExpression(SEXP quantity,SEXP FSE){
    return wrap(
		ExpressionFactory::getInstance()->
		generateXPression(
				  as<string>(quantity),
				  as<string>(FSE) 
				  )
		);
  }



};
