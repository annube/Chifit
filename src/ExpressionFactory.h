


#ifndef EXPRESSION_FACTORY
#define EXPRESSION_FACTORY

#include <map>
#include <string>

using namespace std;

#include <ginac/ginac.h>

using namespace GiNaC;

#include <Rcpp.h>

using namespace Rcpp;


#include "symbols.h"


namespace chifit {


  typedef map<int,ex> ExFromIdMap;
  typedef map<int,ParameterMap> PMapFromIdMap;

  class ExpressionFactory {

  private:

    /**
     * instance to ExpressionFactory
     */
    static ExpressionFactory *instance;

    /**
     * a map storing all expressions by id
     */
    ExFromIdMap XMap;
    PMapFromIdMap PMMap;

    /**
     * the next id to generate
     */
    int nextId;


    /**
     * return the next id to generate
     * use only this function to access nextId
     * accept the constructor (where nextId is only initialized)
     */
    int getNextId(){ return nextId++; }


  public:
    /**
     * return an expression previously generated
     */
    const ex & getExpression(int id){ return XMap[id]; }

    /**
     * return a parameter map for a prev. generated expression from an id
     */
    const ParameterMap & getParameterMap(int id){ return PMMap[id]; }


  private:
    /**
     * constructor
     */
    ExpressionFactory();

  public:

    /**
     * this is the publicly callable function for getting an
     * instance of the factory and should be called instead
     * of the constructor
     * this ensures that we have only one global instance of
     * this class
     */
    static ExpressionFactory * getInstance();

    /**
     * generate expressions
     */
    int generateXPression(const string &quant,const string &FSE);


  };


  /**
   * Rcpp wrapper function for generating expressions
   */
  RcppExport SEXP GenerateChifitExpression(SEXP quantity,SEXP FSE);

};


#endif
