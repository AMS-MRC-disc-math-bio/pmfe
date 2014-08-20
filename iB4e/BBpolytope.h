// Copyright 2006 Peter Huggins
// Released under the GNU GPL license

/*
    This file is a part of iB4e, which is free software;
    you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/



#include "euclideanvector.h"
#include "faces.h"
#include "stack.h"
#include "config.h"
#include <vector>
#include <set>


mpq_class gcd(mpq_class *myvector, int size);
mpq_class gcd2(mpq_class a, mpq_class b);


class vertexnode {
 public:
  vertexnode *next;
  int item;
};



class BBPolytope {
 public:
  BBPolytope(int dim);
  void Build();
  virtual EuclideanVector* BlackBoxOptimize(EuclideanVector *objective) = 0;
  std::vector<EuclideanVector> vertices;

 private:
  //virtual void printNumber(mpq_class a) = 0;
  void processhorizonridges(EuclideanVector *, Face *, Stack *, Stack *); //recursively adds new faces for all horizon ridges, and also marks visible faces as deleted
  Face * vertexandridge(EuclideanVector *v, Face *r);

  void pushvertexintoincidence(int location, EuclideanVector *vertex);
  void printNormals(Face *myface);
  void populateVertices(Face *myface);
  void populateVertexSet(Face *myface, std::set<EuclideanVector> *vertexset);
  void printIncidences(); 
  bool newdirection(Face *myface);
  bool hash(mpq_class *myvector, Face *myface, int recordvertices);

  mpq_class **hashtable;
  int numvertices;
  Stack stack, stack2, stack3, tobedeletedstack;
  Face *dequeue, *dequeue2, *neighbor;
  int dimension;
  vertexnode **facetvertextable;
  int linealitydim;
  int *ghostvertices;

};
