


#include "symbols.h"

#include <ginac/ginac.h>

using namespace GiNaC;

#include <vector>
#include <map>


namespace chifit {

  symbol B("B_0");
  symbol f("f");
  symbol mu("\\mu_q");
  symbol Lambda1("\\Lambda_1");
  symbol Lambda2("\\Lambda_2");
  symbol Lambda3("\\Lambda_3");
  symbol Lambda4("\\Lambda_4");
  symbol Xi3("\\Xi_3");
  symbol L("L");
  symbol ZP("Z_P");
  symbol CMpm("C_{M_{\\pm}}");
  symbol Cf("D_{f}");
  symbol CM0("C_{M_0}");
  symbol c2("c_2");



   std::vector<symbol> allSymbolsOrdered;
   std::vector<double> allSymbolsOrderedDimensions;

  void initAllParamsOrdered(){
    if( allSymbolsOrdered.size() == 0 ){
      allSymbolsOrdered.push_back(B);
      allSymbolsOrderedDimensions.push_back(1);

      allSymbolsOrdered.push_back(f);
      allSymbolsOrderedDimensions.push_back(1);


      allSymbolsOrdered.push_back(c2);
      allSymbolsOrderedDimensions.push_back(4);


      allSymbolsOrdered.push_back(Lambda1);
      allSymbolsOrderedDimensions.push_back(1);

      allSymbolsOrdered.push_back(Lambda2);
      allSymbolsOrderedDimensions.push_back(1);


      allSymbolsOrdered.push_back(Lambda3);
      allSymbolsOrderedDimensions.push_back(1);

      allSymbolsOrdered.push_back(Lambda4);
      allSymbolsOrderedDimensions.push_back(1);

      allSymbolsOrdered.push_back(Xi3);
      allSymbolsOrderedDimensions.push_back(1);

      allSymbolsOrdered.push_back(CMpm);
      allSymbolsOrderedDimensions.push_back(2);

      allSymbolsOrdered.push_back(Cf);
      allSymbolsOrderedDimensions.push_back(2);

      allSymbolsOrdered.push_back(CM0);
      allSymbolsOrderedDimensions.push_back(2);

    }
  }




  /// 

  bool ParameterMap::parUsed(GiNaC::symbol s){
      IndexMapIt itFound = indexMap.find( s );

      if( itFound != indexMap.end() )
	return map[ itFound->second ].second;
      else
	std::cerr << "Error in ParameterMap::parUsed : symbol "<< s <<" not in global symbol table !!" << std::endl;
    }




  ParameterMap::ParameterMap(){
      initAllParamsOrdered();

      for( int i = 0 ; i< allSymbolsOrdered.size(); i++ ){
	map.push_back(std::pair<GiNaC::symbol,bool>(allSymbolsOrdered[i],false));
	indexMap.insert(std::pair<GiNaC::symbol,int>( allSymbolsOrdered[i] , i ) ) ;
      }
    }


    /**
     * mark a parameter as used
     */
  void ParameterMap::add(GiNaC::symbol s){

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

  void ParameterMap::print(){
    for( int i = 0 ; i < map.size() ; i++)
      std::cout << map[i].first << " " << map[i].second << std::endl;
  }






};
