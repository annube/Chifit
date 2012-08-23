#include <ginac/ginac.h>

#include <iostream>
#include <string>
#include <Rcpp.h>

#include <vector>

using namespace std;
using namespace GiNaC;


#include "common.h"


symbol CMpm("C_{M^\\pm}");


ex Mpm_sq = 2*B0*amu_q/a; 
ex M0_sq = 2*B0*amu_q/a; 


ex mpisq = Mpm_sq * ( 1 
		      + M0_sq / 2 / pow( 4 * Pi * f , 2 )  * log( M0_sq / pow(Lambda3,2) ) 
		      + CMpm * pow( a , 2) )  ;


typedef map< string  , vector<double> > AargsType;


RcppExport SEXP mpi_mq_ob(SEXP x, SEXP par,SEXP aargs) {

  /* convert SEXP's to useful Rcpp objects */
   Rcpp::NumericVector vpar(par);
   Rcpp::NumericVector vx(x);
   Rcpp::NumericVector vres(vx.size());
   Rcpp::List aargsMap ( aargs);
   Rcpp::NumericVector latSpac = aargsMap["a"] ;

   /* create a Symbol to input parameter mapping */
   exmap map;
   map[B0] = vpar[0];
   map[f] = vpar[1];
   map[Lambda3] = vpar[2];
   map[CMpm] = vpar[3];

   /* apply the map to the expression */
   ex mpisq_num_eval = evalf( mpisq.subs( map ) );


   /* apply for each value of the vx array and a to the expression and store result in vres*/
   for(int i = 0; i< vx.size() ; i++){
     vres[i] = ex_to<numeric>( 
			      mpisq_num_eval.subs( lst( amu_q == vx[i] ,
							a == latSpac[i] )) 
			      ) .
                      to_double();
   }


   return vres;


}



template <typename T> class setSize 
{

private:
  int size;

public:
  setSize(int n):size(n){}

  void operator()( T &obj ){ obj.reserve(size); }

};


typedef vector<symbol*>::iterator exmapIt;


RcppExport SEXP dmpi_mq_ob(SEXP x, SEXP par,SEXP aargs) {

   Rcpp::NumericVector vpar(par);
   Rcpp::NumericVector vx(x);
   Rcpp::List aargsMap ( aargs);
   Rcpp::NumericVector latSpac = aargsMap["a"] ;


   Rcpp::NumericMatrix gradRes(vx.size(),vpar.size());


   /* create a Symbol to input parameter mapping */
   exmap map;
   map[B0] = vpar[0];
   map[f] = vpar[1];
   map[Lambda3] = vpar[2];
   map[CMpm] = vpar[3];

   vector<symbol*> parvec;
   parvec.push_back(&B0);
   parvec.push_back(&f);
   parvec.push_back(&Lambda3);
   parvec.push_back(&CMpm);



   int jc=0;
   for( exmapIt j=parvec.begin();j!=parvec.end();j++,jc++){

     ex dmpisq = mpisq.diff( *(*j) , 1 );
     ex dmpisq_eval_pars = evalf( dmpisq.subs( map ) );

     //     cout << dmpisq_eval_pars << endl;

     for(int i = 0; i< vx.size() ; i++){

       //       cout << jc << endl;
       gradRes(i,jc) =
	 ex_to<numeric>( dmpisq_eval_pars.subs( lst(
					  amu_q == vx[i],
					  a == latSpac[i]
					  )
				      )	 ).to_double() ;
	      

     }

   }




   return Rcpp::wrap( gradRes );


}
