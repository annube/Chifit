

#include <ginac/ginac.h>

using namespace GiNaC;

#include "../../src/tmFS.h"

int main(int argc,char **argv){

  chifit::ParameterMap pm;
  ex X = get_tm_FSE(pm);

  return 0 ;
}
