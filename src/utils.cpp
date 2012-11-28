


#include "utils.h"

int factorial(int n) {
  if( n == 1 || n == 0 ) return 1;
  int i = n-1,res=n;
  while ( i > 0 )
    res *=i--;
  return res;
}
