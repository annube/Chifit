

#include <vector>
#include <iostream>
#include <iomanip>

using namespace std;


int main(int argc,char** argv){

  vector<double> v(11);

  v[10] = 243;
  v[0] = 3;

  for( int i = 0 ; i < v.size();i++)
    cout << v[i] << endl;

  return 0;
}
