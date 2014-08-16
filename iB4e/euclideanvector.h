


#include <stdlib.h>
#include <iostream>
#include <stdlib.h>

#ifndef EUCLIDEANVECTOR_H
#define EUCLIDEANVECTOR_H

#include "config.h"

class EuclideanVector {

 public:
  int dimension;
  mpq_class *data;

  int id;

  EuclideanVector() {dimension=-1;};
  //EuclideanVector(int d) { dimension = d; data = new mpq_class[d]; for(int i=0; i<d; i++) data[i] = 0;};
  EuclideanVector(int d) { Constructor(d);};
  EuclideanVector(EuclideanVector *v) {
    dimension = v->dimension;
    data = new mpq_class[dimension];
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
  void timesequals(mpq_class c);
  
};

bool computesubfactors(EuclideanVector *orthobasis, int ambientdimension);
mpq_class dotproduct(EuclideanVector *a, EuclideanVector *b);

mpq_class fourbyfourdet(mpq_class a, mpq_class b, mpq_class c, mpq_class d, mpq_class e, mpq_class f, mpq_class g, mpq_class h, mpq_class i, mpq_class j, mpq_class k, mpq_class l, mpq_class m, mpq_class n, mpq_class o, mpq_class p);

//mpq_class threebythreedet(mpq_class a, mpq_class b, mpq_class c, mpq_class d, mpq_class e, mpq_class f, mpq_class g, mpq_class h, mpq_class i); 


#endif

