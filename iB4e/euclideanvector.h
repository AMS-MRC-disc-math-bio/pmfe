


#include <stdlib.h>
#include <iostream>
#include <stdlib.h>

#ifndef EUCLIDEANVECTOR_H
#define EUCLIDEANVECTOR_H

#include "config.h"

class EuclideanVector {

 public:
  int dimension;
  NUMBER *data;

  int id;

  EuclideanVector() {dimension=-1;};
  //EuclideanVector(int d) { dimension = d; data = new NUMBER[d]; for(int i=0; i<d; i++) data[i] = 0;};
  EuclideanVector(int d) { Constructor(d);};
  EuclideanVector(EuclideanVector *v) {
    dimension = v->dimension;
    data = new NUMBER[dimension];
    for(int i =0; i < dimension; i++)
      data[i] = v->data[i];
  };

  void deletedata() { 
     delete[] data;
    return; 
  };

  ~EuclideanVector(){return;};
  
  void Print();
  void Constructor(int d);
  void Negation();
  void plusequals(EuclideanVector w);
  void minusequals(EuclideanVector w);
  void timesequals(NUMBER c);
  
};

bool computesubfactors(EuclideanVector *orthobasis, int ambientdimension);
NUMBER dotproduct(EuclideanVector *a, EuclideanVector *b);

NUMBER fourbyfourdet(NUMBER a, NUMBER b, NUMBER c, NUMBER d, NUMBER e, NUMBER f, NUMBER g, NUMBER h, NUMBER i, NUMBER j, NUMBER k, NUMBER l, NUMBER m, NUMBER n, NUMBER o, NUMBER p);

//NUMBER threebythreedet(NUMBER a, NUMBER b, NUMBER c, NUMBER d, NUMBER e, NUMBER f, NUMBER g, NUMBER h, NUMBER i); 


#endif

