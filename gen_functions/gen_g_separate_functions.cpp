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

  ofstream ofs_h("tmFS.R.g.h");

  ofs << "#include <complex>" << endl;

  ofs_h << "#include <complex>" << endl;

  ofs << "#include <iostream>" << endl;
  ofs << "#include <iomanip>" << endl;

  ofs << "namespace " << namespaceGen << " { " << endl; 

  ofs_h << "namespace " << namespaceGen << " { " << endl; 


  for(int d = 0 ; d <=maxDeri ;d ++) {
    ofs << " std::complex<double> "<< functionName << "_d" << d <<"(std::complex<double> x) { " << endl;
    ofs << "  return( " << endl;
    ofs << csrc_double << g.diff(x,d) << endl;
    ofs << " ) ; " << endl;

    ofs << "}" << endl;

    ofs_h << " std::complex<double> "<< functionName << "_d" << d <<"(std::complex<double> x) ; " << endl;

  }

  ofs << "};" << endl;
  ofs_h << "};" << endl;


  ofs.close();
  ofs_h.close();

}
