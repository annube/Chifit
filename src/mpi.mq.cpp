#include <ginac/ginac.h>

#include <iostream>
#include <string>
#include <Rcpp.h>

#include <vector>

using namespace std;

using namespace GiNaC;





// mpi.mq <- function(x,par,aargs){
//   two.B0 <- par[1]
//   f0 <- par[2]
//   Lambda3 <- par[3]
//   a <- aargs$a
//   Mpm.sq = two.B0*x/a
//   M0.sq =  two.B0*x/a
//   return(  Mpm.sq * ( 1 + M0.sq/2/(4*pi*f0)^2 * log( M0.sq /  Lambda3^2 ) + par[4] * a^2 ) )
  
// }



symbol B0("B_0");
symbol f("f");
symbol a("a");
symbol amu_q("\\mu_q");
symbol Lambda3("\\lambda_3");
symbol CMpm("C_{M^\\pm}");


ex Mpm_sq = 2*B0*amu_q/a; 
ex M0_sq = 2*B0*amu_q/a; 


ex mpisq = Mpm_sq * ( 1 
		      + M0_sq / 2 / pow( 4 * Pi * f , 2 )  * log( M0_sq / Lambda3 ) 
		      )  ;


typedef map< string  , vector<double> > AargsType;


RcppExport SEXP mpi_mq(SEXP x, SEXP par,SEXP aargs) {

   Rcpp::NumericVector vpar(par);
   Rcpp::NumericVector vx(x);
   Rcpp::NumericVector vres(vx.size());

   Rcpp::List aargsMap ( aargs);

   Rcpp::NumericVector latSpac = aargsMap["a"] ;


   exmap map;
   map[B0] = vpar[0];
   map[f] = vpar[1];
   map[Lambda3] = vpar[2];



   ex mpisq_num_eval = evalf( mpisq.subs( map ) );


   for(int i = 0; i< vx.size() ; i++){

     vres[i] = ex_to<numeric>( 
			      mpisq_num_eval.subs( mu_q == vx[i] ,
						   a == latSpac[i] ) 
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



RcppExport SEXP dmpi_mq(SEXP x, SEXP par,SEXP aargs) {

   Rcpp::NumericVector vpar(par);
   Rcpp::NumericVector vx(x);
   Rcpp::NumericVector vres(vx.size());


   Rcpp::NumericMatrix gradRes(vx.size(),vpar.size());


   symbol B0("B_0");
   symbol f("f");
   symbol mu_q("\\mu_q");

   vector<ex> parExs(2);

   ex mpisq = 2*B0*mu_q * ( 1 + 1/2);

   ex dmpisq_dB0 = mpisq.diff(B0,1);

   ex dmpisq_dB0_num_eval = evalf( dmpisq_dB0.subs( B0 == vpar[0] ) );



   for( int j = 0 ; j < vpar.size() ; j++ ){

     for(int i = 0; i< vx.size() ; i++){

       gradRes(i,j) =
	 ex_to<numeric>( dmpisq_dB0_num_eval.subs( mu_q == vx[i] ) ).to_double() ;
	      

     }
   }




   return Rcpp::wrap( gradRes );


}
