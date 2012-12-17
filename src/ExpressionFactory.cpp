
#include <iostream>

using namespace std;


#include <ginac/ginac.h>

using namespace GiNaC;

#include <Rcpp.h>

using namespace Rcpp;

#include "symbols.h"
#include "mpi.mq.ob.h"
#include "tmFS.h"
#include "mpi.mq.ob.fse.h"
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

    int id = getNextId();

    if( quant == "mpi.mq.ob" ){
      ParameterMap pm;

      XMap[id] = get_M_pm_sq_Xpression(pm,false);

      if( FSE == "tmFS" )
	XMap[id] = XMap[id] * get_tm_FSE_Mpm_sq( pm );
      else if(FSE == "ob" )
	XMap[id] = XMap[id] + get_M_pm_sq_FSE_ob_Xpression( pm );



      PMMap[id] = pm;
    }

    return id;
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
