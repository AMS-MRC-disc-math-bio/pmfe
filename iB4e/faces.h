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




#ifndef FACES_H
#define FACES_H


#include "euclideanvector.h"
#include "config.h"
#include <iostream>
#include <stdlib.h>

#define MAXFACETS   1000000
#define MAXVERTICES 5000000

using namespace std;

class Face {

 public:
  int dimension;
  int ambientdimension;
  Face **facets;
  int maxfacets;
  int numfacets;
  EuclideanVector **normalvectors;
  mpq_class *rhs;
  Face **incidents;
  EuclideanVector **vertices;
  bool deleted;
  int maxvertices;
  int numvertices;
  bool enqueued;
  //int marker;
  bool tobedeleted;

  Face() { 
                facets = NULL;
                maxfacets = 0;
                numfacets = 0;
                normalvectors = NULL; 
                rhs = NULL; 
                dimension = -2;
                ambientdimension = -2; 
                incidents = NULL; 
                vertices = NULL;
                maxvertices = 0;
                numvertices = 0; 
                deleted = false;
                enqueued = false;
                tobedeleted = false;
			   // marker = 0;	
              };

  void deleteme();
  void InitializeMaxFacets(int f); // Run when a new face is created
  void InitializeDimensions(int dim, int ambientdim);  // Run when a new face is created
  void InitializeMaxVertices(int v); // Run when a new face is created
  void AddLIVertexToLowDim(EuclideanVector *v);  // Assumes that vertex is not in the affine span of the rest! 
  Face * AddFacetWithoutVertex(int i);
  void SortVertices();
  void DeleteFacet(Face *f);
  void AddFacet(Face *f);
  void ComputeNormalVectorAwayFromPoint(EuclideanVector *p);  
  void AddVertex(EuclideanVector *v);
  void ComputeRHS();
  void PrintVertices();
  void Print();

};

#endif
