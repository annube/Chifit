
#include <iostream>

using namespace std;

#define MATHLIB_STANDALONE 

//extern "C" {
#include "Rmath.h"
//}

int main(int argc, char ** argv){
  cout << bessel_k( 5. , 1. , 1.) << endl;
  return 0;
}
