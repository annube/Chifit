

#include <ginac/ginac.h>

using namespace GiNaC;

#include <iostream>

using namespace std;

DECLARE_FUNCTION_1P(gtilde);

DECLARE_FUNCTION_1P(gtilde_deri);


static ex gtilde_deri_eval(const ex &x){
  return gtilde_deri(x).hold();
}

static ex gtilde_deri_evalf(const ex &x){

  if(is_a<numeric>(x) ){
    return 1;
  }
  else 
    return gtilde_deri(x).hold();
}

static ex gtilde_deri_eval_deri(const ex &x, unsigned diffpar){
  return 1/0;
}

REGISTER_FUNCTION(gtilde_deri,eval_func(gtilde_deri_eval).
		  evalf_func(gtilde_deri_evalf).
		  derivative_func(gtilde_deri_eval_deri)
		  );


static ex gtilde_eval(const ex &x){
  return gtilde(x).hold();
}

static ex gtilde_evalf(const ex &x){

  if(is_a<numeric>(x) ){
    return 1;
  }
  else 
    return gtilde(x).hold();
}

static ex gtilde_eval_deri(const ex &x, unsigned diffpar){
  return gtilde_deri(x);
}


REGISTER_FUNCTION(gtilde,eval_func(gtilde_eval).
		  evalf_func(gtilde_evalf).
		  derivative_func(gtilde_eval_deri)
		  );

int main(int argc,char** argv){

  symbol x("x");

  ex formula=x*gtilde(x);

  cout << formula << endl;

  ex dformula_dx=formula.diff(x,1);

  cout << dformula_dx << endl;



  return 0;
}
