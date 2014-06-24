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



#include "faces.h"
#include <iostream>
#include <stdlib.h>
#include <stdint.h>

//#define DEBUG

using namespace std;


void Face::InitializeMaxFacets(int f) 
{

  maxfacets = f;
  numfacets = 0;
  facets = new Face *[f];
// MAY NEED TO UNCOMMENT!!
  //rhs = new double[f];

  return;
}


void Face::InitializeMaxVertices(int v) 
{

  maxvertices = v;
  numvertices = 0;
  vertices = new EuclideanVector *[v];

  return;
}

void Face::InitializeDimensions(int dim, int ambientdim)
{

  if (dim > ambientdim) {
    cerr << "dimension greater than ambient dimension!";
    return;
  }

  if (dim < -1 || ambientdim < 1) {
    cerr << "dimension or ambient dimension out of range!";
    return;
  }

  dimension = dim;
  ambientdimension = ambientdim;

  if(ambientdim > dim) {
    if(dim > -1) {
      normalvectors = new EuclideanVector *[ambientdim - dim];
      rhs = new NUMBER[ambientdim - dim];
      for(int i = 0; i < ambientdim - dim; i++)
        normalvectors[i] = NULL;
    }
    else { 
      normalvectors = new EuclideanVector *[ambientdim];
      rhs = new NUMBER[ambientdim];
      for(int i = 0; i < ambientdim; i++)
        normalvectors[i] = NULL;
    }
  }

  for(int i = 0; (i < ambientdim - dim) && (i < ambientdim) ; i++)  {
    normalvectors[i] = new EuclideanVector(ambientdim);
  }

  return;

}

void Face::deleteme()
{

  if(normalvectors != NULL)
  for(int i = 0; (i < ambientdimension - dimension) && (i < ambientdimension) ; i++) {
    
    if(normalvectors[i] != NULL) {
      normalvectors[i]->deletedata();
      delete normalvectors[i];
    }
  }
   
  if(rhs != NULL)
    delete[] rhs;
 
  if(incidents != NULL)
    delete[] incidents;
 
  if(normalvectors != NULL)
    delete[] normalvectors;
  if(vertices != NULL)
    delete[] vertices; 
  if(facets != NULL)
    delete[] facets;


  return;
} 


void Face::Print()
{

  #ifdef DEBUG2
  cout << "Face::Print() called...\n";
  #endif
  cout << "Address of this Face:  " << this << "\n";
  cout << "Vertices of " << this << ":\n";
  PrintVertices();
  
  cout << "Ridges and neighbors of " << this << ":\n";
  for(int i = 0; i < numfacets; i++) {
    cout << "Ridge " << i << ": " << facets[i];
    cout << "Incidences to ridge " << i << ":\n";
    cout << "   " << facets[i]->incidents[0] << "   " << facets[i]->incidents[1];    cout << "\nVertices of ridge " << i << ":\n";
    facets[i]->PrintVertices(); 
  }  

  return;
}


void Face::PrintVertices()
{
  #ifdef DEBUG2
  cout << "\nFace::PrintVertices() called...Vertices are: ";
  #endif

  for(int i = 0; i < numvertices; i++)
    vertices[i]->Print();

}

void Face::AddLIVertexToLowDim(EuclideanVector *v)
{
  dimension++;
  //vertices[numvertices] = new EuclideanVector(ambientdimension);
  vertices[numvertices] = v;

  //for(int i = 0; i < ambientdimension; i++)
  //  vertices[numvertices]->data[i] = v->data[i];
  
  numvertices++;

  //cout << "\nambient dimension: " << ambientdimension << "\n";

  if(numvertices == 1) {
    for(int i = 0; i < ambientdimension; i++)
    for(int j = 0; j < ambientdimension; j++)  {
        if( i == j)
          normalvectors[i]->data[j] = 1;
        else
          normalvectors[i]->data[j] = 0;
    }
  }

  if((numvertices) > 1 && (numvertices <= ambientdimension)) {
    //EuclideanVector orthobasis[ambientdimension];
    EuclideanVector *orthobasis = new EuclideanVector[ambientdimension];

    for(int i = 0; i < ambientdimension; i++)
      orthobasis[i].Constructor(ambientdimension);

    for(int i = 1; i < numvertices; i++) {
      orthobasis[i-1].plusequals(*(vertices[0]));
      orthobasis[i-1].timesequals(-1);
      orthobasis[i-1].plusequals(*(vertices[i]));
      #ifdef DEBUG
      cout << "Add2LI:  Orthobasis vector :  ";
      orthobasis[i-1].Print();
      #endif 
    }

    bool foundnullvector = false;

    while(!foundnullvector) {

      #ifdef DEBUG
        cout << "Trying to find a null vector in lowdim\n";
        cout << "with " << numvertices << "  vertices\n";
      #endif
    
      for(int i = numvertices - 1; i < ambientdimension; i++) {
        for(int j = 0; j < ambientdimension; j++)
          orthobasis[i].data[j] = rand()  % 103;
        //orthobasis[i].Print(); 
      } 

      foundnullvector = computesubfactors(orthobasis, ambientdimension);

      #ifdef DEBUG
      if(!foundnullvector) {
        cout << "whoops we didn't find a nullvector...\n";
        PrintVertices();
      }
      else {
        cout << "Found nullvector:  ";
        orthobasis[ambientdimension - 1].Print();
      } 

      #endif

    }

    for(int i = 0 ; i < ambientdimension; i++) 
      normalvectors[0]->data[i] = orthobasis[ambientdimension - 1].data[i];

  }

  if(numvertices == ambientdimension + 1) {  // we have a simplex

    #ifdef DEBUG
    cout << "\n AddLowDim:  Houston we have a simplex";
    #endif  

    //exit(0);

    EuclideanVector centroid(ambientdimension);
    
    for(int i = 0; i < ambientdimension + 1; i++) 	  
      centroid.plusequals(*(vertices[i]));
    
    
    for(int i = 0; i < ambientdimension+1; i++) {
      numfacets++;
      facets[i] = new Face;
      facets[i]->InitializeDimensions(-1, ambientdimension);
      facets[i]->InitializeMaxFacets(ambientdimension);
      facets[i]->InitializeMaxVertices(ambientdimension);
      
      int differentj;
     
 
      for(int j = 0; j < ambientdimension+1; j++)
        if(j != i) {
	  facets[i]->AddLIVertexToLowDim(vertices[j]);
	  differentj = j;
	}


      facets[i]->rhs[0] = dotproduct(facets[i]->normalvectors[0], vertices[differentj]);

      if((ambientdimension+1)*facets[i]->rhs[0] < dotproduct((facets[i]->normalvectors[0]), &centroid)) {
	   (facets[i]->normalvectors[0])->Negation();
           facets[i]->rhs[0] *= -1;
      }

      #ifdef DEBUG
      cout << "\nNormal vector for facet " << i << " is ";
      facets[i]->normalvectors[0]->Print();
      cout << "with RHS = " << facets[i]->rhs[0];
      #endif

      facets[i]->SortVertices();
      
    }  //end for

    
    Face *newridge;

    for(int i = 0; i < numfacets; i++)
    for(int j = 0; j < i; j++) {
      #ifdef DEBUG
      cout <<  "\nWe're creating the first ridges!";
      #endif
      newridge = new Face;

      //debug
      //return;

      newridge->InitializeDimensions(dimension-2, ambientdimension);
      newridge->InitializeMaxVertices(dimension-1);

      for(int k = 0; k < numvertices; k++)
      if(k != i)
      if(k != j)
        newridge->AddVertex(vertices[k]);

      facets[i]->AddFacet(newridge);
      facets[j]->AddFacet(newridge);
      
      newridge->incidents = new Face *[2];

      newridge->incidents[0] = facets[i];
      newridge->incidents[1] = facets[j];
    }

  }

  return;

}  

void Face::DeleteFacet(Face *f)
{
  //probably this is leading to a slow memory leak

  int i;

  for( i = 0; facets[i] != f; i++)
    ; // NULL FOR

  facets[i]->deleteme();

  facets[i] = facets[numfacets - 1];
  numfacets--;

  return;

}

void Face::AddFacet(Face *f)
{
  facets[numfacets] = f;
  numfacets++;

  return;
}

void Face::AddVertex(EuclideanVector *v)
{
  vertices[numvertices] = v;
  numvertices++;

  return;

}

Face *  Face::AddFacetWithoutVertex(int k)
{
  facets[numfacets] = new Face;
  
  //facets[numfacets]->InitializeMaxFacets(dimension);
  facets[numfacets]->InitializeDimensions(dimension-1, ambientdimension);
  facets[numfacets]->InitializeMaxVertices(dimension);

  for(int i = 0; i < numvertices; i++)
    if(i != k)
      facets[numfacets]->AddVertex(vertices[i]);

  numfacets++;

  return facets[numfacets - 1];
}

void Face::ComputeNormalVectorAwayFromPoint(EuclideanVector *p) 
{

  //EuclideanVector orthobasis[ambientdimension];
  EuclideanVector *orthobasis = new EuclideanVector[ambientdimension];

  for(int i = 0; i < ambientdimension; i++)
    orthobasis[i].Constructor(ambientdimension);
    
  for(int i = 1; i < numvertices; i++) {
  
      orthobasis[i-1].plusequals(*(vertices[0]));
      orthobasis[i-1].timesequals(-1);
      orthobasis[i-1].plusequals(*(vertices[i]));

  }
    
  bool foundnullvector = false;

  foundnullvector = computesubfactors(orthobasis, ambientdimension);

  if(!foundnullvector) {
    #ifdef DEBUG
    PrintVertices();
    #endif
    cerr << "Facet not full-dimensional!!!\n";
    exit(0);
  }

  for(int i = 0; i < ambientdimension; i++)
    normalvectors[0]->data[i] = orthobasis[ambientdimension - 1].data[i];

  rhs = new NUMBER[1];
  *rhs = dotproduct(normalvectors[0], vertices[0]);

  if((ambientdimension+1)*(*rhs) < dotproduct(p, normalvectors[0])) {
    *rhs = -1 * (*rhs);
    normalvectors[0]->Negation();
  }

  return;
}

void Face::SortVertices()
{
  EuclideanVector *swap;

  uintptr_t compare1, compare2;

  for(int i = 0; i < numvertices; i++)
  for(int j = 0; j < numvertices; j++)
  {
    compare1 = (uintptr_t)((void *)(vertices[i]));
    compare2 = (uintptr_t)((void *)(vertices[j]));
    if(compare1 < compare2) {
      swap = vertices[i];
      vertices[i] = vertices[j];
      vertices[j] = swap;
    }

  }

  return;

}

void Face::ComputeRHS()
{

  for(int i = 0; i < ambientdimension - dimension; i++) {
    rhs[i] = dotproduct(vertices[numvertices-1], normalvectors[i]);

    for(int j = 0; j < numvertices-1; j++)
      if(rhs[i] < dotproduct(vertices[j], normalvectors[i]))
        rhs[i] = dotproduct(vertices[j], normalvectors[i]);

  }

}
