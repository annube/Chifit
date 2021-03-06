
#include <iostream>

using namespace std;


#include <ginac/ginac.h>

using namespace GiNaC;

#include <Rcpp.h>

using namespace Rcpp;

#include "symbols.h"
#include "mpi.mq.ob.h"
#include "mpi.mq.ob.fse.h"
#include "tmFS.h"
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

    ParameterMap pm;


    if( quant == "mpi.mq.ob" ){

      XMap[id] = get_M_pm_sq_ob_Xpression(pm,false);

      if(FSE == "ob" )
	XMap[id] = XMap[id] + get_M_pm_sq_ob_FSE_Xpression( pm );
      else if( FSE == "tmFS" )
	XMap[id] = XMap[id] * pow( 1 + get_tm_FSE_Rpm( pm ) , 2 );

      PMMap[id] = pm;

    } else if( quant == "fpi.mq.ob" ){

      XMap[id] = get_f_pm_ob_Xpression(pm,false);

      if(FSE == "ob" )
	XMap[id] = XMap[id] + get_f_pm_ob_FSE_Xpression( pm );
      else  if(FSE == "tmFS" )
	XMap[id] = XMap[id] * ( 1 +  get_tm_FSE_Rfpm( pm  ) );


      PMMap[id] = pm;
    } else if( quant == "mpi_0.mq.ob" ){

      XMap[id] = get_M_0_sq_ob_Xpression(pm,false);

      if(FSE == "ob" )
	XMap[id] = XMap[id] + get_M_0_sq_ob_FSE_Xpression( pm );
      else  if(FSE == "tmFS" )
	XMap[id] = XMap[id] * pow( 1 +  get_tm_FSE_RM0( pm ) , 2 );


      PMMap[id] = pm;
    } else  if( quant == "ob.FSE.mpm" ) {
      XMap[id] =  get_M_pm_sq_ob_FSE_Xpression( pm ) ;
      PMMap[id] = pm;
    } else  if( quant == "ob.FSE.fpm" ) {
      XMap[id] = get_f_pm_ob_FSE_Xpression( pm );
      PMMap[id] = pm;
    } else  if( quant == "tmFS.mpm" ) {
      XMap[id] =  1. + 2. * get_tm_FSE_Rpm( pm );
      PMMap[id] = pm;
    } else  if( quant == "tmFS.fpm" ) {
      XMap[id] = 1. +  get_tm_FSE_Rfpm( pm );
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
