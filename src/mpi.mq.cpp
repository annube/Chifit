#include <ginac/ginac.h>

#include <iostream>
#include <string>
#include <Rcpp.h>

#include <vector>

using namespace std;

using namespace GiNaC;



#include "common.h"


static ex chi_mu = 2*B0*amu_q/a; 


static ex ampisq = pow(a , 2) * chi_mu * ( 1 
				    + pow(a,  2) * chi_mu / pow( 4 * Pi * a * f , 2 )  * log( pow(a,2)*chi_mu / pow( a * Lambda3,2) ) 
			       )  ;


typedef map< string  , vector<double> > AargsType;


RcppExport SEXP mpi_mq(SEXP x, SEXP par,SEXP aargs) {

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

   /* apply the map to the expression */
   ex ampisq_num_eval = evalf( ampisq.subs( map ) );


   /* apply for each value of the vx array and a to the expression and store result in vres*/
   for(int i = 0; i< vx.size() ; i++){
     vres[i] = ex_to<numeric>( 
			      ampisq_num_eval.subs( lst( amu_q == vx[i] ,
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


typedef vector<const symbol*>::iterator exmapIt;
// typedef vector<const symbol*>::const_iterator exmapCIt;


RcppExport SEXP dmpi_mq(SEXP x, SEXP par,SEXP aargs) {

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

   vector<const symbol*> parvec;
   parvec.push_back(&B0);
   parvec.push_back(&f);
   parvec.push_back(&Lambda3);



   int jc=0;
   for( exmapIt j=parvec.begin();j!=parvec.end();j++,jc++){

     ex dampisq = ampisq.diff( *(*j) , 1 );
     ex dampisq_eval_pars = evalf( dampisq.subs( map ) );

     //     cout << dmpisq_eval_pars << endl;

     for(int i = 0; i< vx.size() ; i++){

       //       cout << jc << endl;
       gradRes(i,jc) =
	 ex_to<numeric>( dampisq_eval_pars.subs( lst(
					  amu_q == vx[i],
					  a == latSpac[i]
					  )
				      )	 ).to_double() ;
	      

     }

   }




   return Rcpp::wrap( gradRes );


}
