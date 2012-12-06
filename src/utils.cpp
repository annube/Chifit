


#include "utils.h"

#include <iostream>

int factorial(int n) {
  if( n == 1 || n == 0 ) return 1;
  int i = n-1,res=n;
  while ( i > 0 )
    res *=i--;
  return res;
}

int MyExToInt(GiNaC::ex i){
  using namespace std;
  using namespace GiNaC;

  int in=(int) ( ex_to<numeric>(i).to_double() ) ;

  double i_fract = ex_to<numeric>(i).to_double() - (double)in;
  if( i_fract != 0.0 ) cerr << "Error i is not integer in MyExToInt" << endl;
  
  return in;
    
}
