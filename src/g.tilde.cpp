#include <iostream>
#include <iomanip>
#include <tr1/cmath>
#include <Rcpp.h>
#include <math.h>
#include <vector>

// #include <boost/mpl.h>
// #include <boost/mpl/eval_if.hpp>

using namespace std;
using namespace tr1;

template <typename T>
T mysin( const T &x){
  return (T)sin((double)x);
}

template <int n>
double mybessel( const double &x){
  return cyl_bessel_k(n,x);
}

typedef double (*applyfn)(const double &) ;

applyfn besseli[]={mybessel<0>,mybessel<1>,mybessel<2>};



int partitions[]={6.,12.,8.,6.,24.,24.,0.,12.,30.,24.,24.,8.,24.,48.,0.,6.,48.,36.,24.,24.};
int numPartitions = sizeof(partitions)/sizeof(int);

double g_tilde_scalar(double x){
  double result=0;
  for( int  i=0;i<numPartitions;i++) {
    double y = x * sqrt(i+1);
    result += 4 * partitions[i]*cyl_bessel_k(1,y)/y ;
  }
  return result;
}

RcppExport SEXP g_tilde(SEXP l,SEXP n){


    using namespace Rcpp ;
    
    NumericVector xl(l);

    NumericVector res(xl.size());
    int index=as<int>(n);

    if(index < 0 | index >2 ) {
      cerr << "Error index out of range" << endl;
      return res;
    }

    res = sapply(xl,g_tilde_scalar);
    return res;
}


// g1 <- function(x) {

//   weights <- c(6.,12.,8.,6.,24.,24.,0.,12.,30.,24.,24.,8.,24.,48.,0.,6.,48.,36.,24.,24.)
//   ex <- c(1:20)
//   res <- x
//   for( i in 1:length(x)) {
//     sex <- x[i]*sqrt(ex)
//     res[i] <- sum(4*weights*besselK(sex, 1)/(sex))
//   }
//   return(res)
// }
