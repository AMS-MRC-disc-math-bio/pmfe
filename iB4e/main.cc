#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <string.h>
#include <vector>
#include <algorithm>

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
  general_Polytope(int d, const char *call):BBPolytope(d){strcpy(optimizer_caller, call);};
  void BlackBoxOptimize(EuclideanVector *objective, EuclideanVector *solution);
  char optimizer_caller[500];
};


int main(const int argc, const char * argv[]) 
{
  int dim;
  sscanf(argv[1],"%d",&dim);

  general_Polytope polytope(dim, argv[2]); // second arg is the call string
  polytope.Build();

  return 0;
}


void general_Polytope::BlackBoxOptimize(EuclideanVector *objective, EuclideanVector *solution) 
{

  char optimizer_call[1000];

  ofstream out("./objective");    

  for(int i = 0; i < dimension; i++) {
    #ifdef GMP
      out << objective->data[i].get_str(10) << " ";
    #else
      out << objective->data[i] << " ";
    #endif
  }
  out.close();
  
  sprintf(optimizer_call, "cat ./objective | %s  > ./solution ", optimizer_caller); 
  system(optimizer_call);

  ifstream in("./solution");
  for(int i = 0; i < dimension; i++)
      in >> solution->data[i]; // does this work with GMP ????
  in.close();

  return;

}
