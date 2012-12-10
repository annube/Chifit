

// link against tmFS.R.g.o
#include <complex>
#include <iostream>

using namespace std;

namespace chifit {
complex<double> R_g(complex<double> x , int deri=0);
}

int main(int argc,char ** argv){

  cout << chifit::R_g( complex<double>(1,0),0 ) << endl;
  cout << chifit::R_g( complex<double>(1,0),1 ) << endl;
  cout << chifit::R_g( complex<double>(1,0),2 ) << endl;
  cout << chifit::R_g( complex<double>(1,0),3 ) << endl;

  return 0;
}
