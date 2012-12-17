


#ifndef EXPRESSION_FACTORY
#define EXPRESSION_FACTORY

#include <map>
#include <string>

using namespace std;

#include <ginac/ginac.h>

using namespace GiNaC;

#include <Rcpp.h>

using namespace Rcpp;


namespace chifit {


  typedef map<int,ex> ExFromIdMap;

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


    /**
     * constructor
     */
    ExpressionFactory();

  public:

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
