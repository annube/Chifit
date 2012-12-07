

#ifndef SYMBOLS_H
#define SYMBOLS_H


#include <ginac/ginac.h>

#include <map>
#include <vector>
#include <iostream>
#include <iomanip>

namespace chifit {

  extern  GiNaC::symbol B;
  extern  GiNaC::symbol f;
  extern  GiNaC::symbol mu;
  extern  GiNaC::symbol Lambda1;
  extern  GiNaC::symbol Lambda2;
  extern  GiNaC::symbol Lambda3;
  extern  GiNaC::symbol Lambda4;
  extern  GiNaC::symbol Xi3;
  extern  GiNaC::symbol L;
  extern  GiNaC::symbol ZP;
  extern  GiNaC::symbol CMpm;
  extern  GiNaC::symbol Cf;
  extern  GiNaC::symbol CM0;
  extern  GiNaC::symbol c2;


  extern std::vector<GiNaC::symbol> allSymbolsOrdered;
  extern std::vector<double> allSymbolsOrderedDimensions;


  void initAllParamsOrdered();


  typedef std::vector<std::pair<GiNaC::symbol,bool> > SymbolBoolVec;
  typedef SymbolBoolVec::iterator SymbolBoolVecIt;


  /**
   * class for creating a map of parameters occuring in certain XPT expressions
   */


  class ParameterMap {

  public:
    typedef std::map<GiNaC::symbol,int,GiNaC::ex_is_less> IndexMap;
    typedef IndexMap::iterator IndexMapIt;

  private:
    SymbolBoolVec map;

    IndexMap indexMap;

  public:

    const SymbolBoolVec & getMap(){ return map; }

    bool parUsed(GiNaC::symbol s){
      IndexMapIt itFound = indexMap.find( s );

      if( itFound != indexMap.end() )
	return map[ itFound->second ].second;
      else
	std::cerr << "Error in ParameterMap::parUsed : symbol "<< s <<" not in global symbol table !!" << std::endl;
    }

  public:

    /**
     * Constructor
     */

    ParameterMap(){
      initAllParamsOrdered();

      for( int i = 0 ; i< allSymbolsOrdered.size(); i++ ){
	map.push_back(std::pair<GiNaC::symbol,bool>(allSymbolsOrdered[i],false));
	indexMap.insert(std::pair<GiNaC::symbol,int>( allSymbolsOrdered[i] , i ) ) ;
      }
    }


    /**
     * mark a parameter as used
     */

    void add(GiNaC::symbol s){

/*       std::cout << "printing index map " << std::endl; */
/*       for(IndexMapIt imi=indexMap.begin() ; imi != indexMap.end() ; imi++) */
/* 	std::cout << imi->first << " " << imi->second << std::endl; */

//      std::cout << "found symbol at index " << indexMap[ s ] << std::endl;
      IndexMapIt itFound = indexMap.find( s );

      if( itFound != indexMap.end() )
	map[ itFound->second ].second = true;
      else
	std::cerr << "Error in ParameterMap::add : symbol "<< s <<" not in global symbol table !!" << std::endl;
    }


    /**
     * print symbol map
     */

    void print(){
      for( int i = 0 ; i < map.size() ; i++)
	std::cout << map[i].first << " " << map[i].second << std::endl;
    }

  };


};

#endif
