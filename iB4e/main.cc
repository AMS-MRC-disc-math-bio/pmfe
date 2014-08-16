#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>

#include "BBpolytope.h"
#include "euclideanvector.h"
#include "faces.h"

#include "config.h"

#ifdef GMP
  #include <gmpxx.h>
#endif

using namespace std;


class general_Polytope : public BBPolytope {

 public:
  general_Polytope(int d):BBPolytope(d){;};
  void BlackBoxOptimize(EuclideanVector *objective, EuclideanVector *solution);
};


int main(const int argc, const char * argv[]) 
{
  int dim;
  sscanf(argv[1],"%d",&dim);

  general_Polytope polytope(dim);
  polytope.Build();

  return 0;
}


void general_Polytope::BlackBoxOptimize(EuclideanVector *objective, EuclideanVector *solution) 
{
  cout << "get_result([";
  for(int i = 0; i < dimension; i++) {
    #ifdef GMP
      cout << objective->data[i].get_str(10) << ",";
    #else
      cout << objective->data[i] << ",";
    #endif
  }
  cout << "])\n";

  string line;
  getline(cin, line);
  cin.clear();

  string buf;
  stringstream ss(line);
  vector<long> tokens;

  while (ss >> buf)
    tokens.push_back(atol(buf.c_str()));
  
  for(int i = 0; i < dimension; i++){
    solution->data[i] = tokens[i];
  }
 
  return;
}
