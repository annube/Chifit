#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include <ginac/ginac.h>

using namespace GiNaC;

int main(int argc,char **argv){

  string namespaceGen("chifit");
  string functionName("R_g");
  const int maxDeri = 2;

  symbol x("x");
  ex sigma = sqrt(1.-4./x);
  ex g=sigma*log( (sigma - 1. ) / (sigma + 1. ) ) + 2.; 


  ofstream ofs("tmFS.R.g.cpp");

  ofs << "#include <complex>" << endl;
  ofs << "#include <iostream>" << endl;
  ofs << "#include <iomanip>" << endl;

  ofs << "namespace " << namespaceGen << " { " << endl; 

  ofs << " std::complex<double> "<< functionName <<"(std::complex<double> x,int deri=0) { " << endl;
  ofs << "switch(deri) { " << endl;

  for(int d = 0 ; d <=maxDeri ;d ++) {
    ofs << " case " << d << ": return( " << endl;
    ofs << csrc_double << g.diff(x,d) << endl;
    ofs << " ) ; " << endl;
    ofs << "break;" << endl;

  }

  ofs << " default: std::cerr << \"Error deri out of range should be <= " << maxDeri << "!!!\" << std::endl; return NAN ; break; " << endl;

  ofs << "}" << endl;

  ofs << "}" << endl;

  ofs << "};" << endl;


  ofs.close();

}
