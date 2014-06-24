
#include <stdlib.h>

#include "stack.h"
#include "faces.h"

//#define DEBUG

void Stack::Push(void *x) 
{
  (blocks[currentblock])[currentpos] = x;
  currentpos++;
  totalsize++;
  if(currentpos == STACKBLOCKSIZE) {
	currentpos = 0;
	currentblock++;
	blocks[currentblock] = new void *[STACKBLOCKSIZE];
  }

  empty = false;

 // ((Face *)(x))->enqueued = true;

  return;
}

void * Stack::Pop() 
{
  if((currentpos == 0) && (currentblock == 0))
	return NULL;
  
  currentpos--;
  totalsize--;
  
  if(currentpos == -1) {
	currentpos = STACKBLOCKSIZE - 1;
	currentblock--;
  }


  if((currentpos == 0) && (currentblock == 0))
	empty = true;


  //((Face *)((blocks[currentblock])[currentpos]))->enqueued = false;
  
  return (blocks[currentblock])[currentpos];

}

void Stack::Sort()
{
   void *swap;

   if(totalblocks > 1) {
     cerr << "We can't sort big stacks yet!";
     exit(1);
   }

   for(int i = 0; i < totalsize; i++)
   for(int j = 0; j < totalsize; j++) {
     if(sortcompare((blocks[0])[i], (blocks[0])[j]) < 0) {
       swap = (blocks[0])[i];
       (blocks[0])[i] = (blocks[0])[j];
       (blocks[0])[j] = swap;
     }
   }
   
  return;
}

int sortcompare(void *a, void *b) 
{

  Face *face1, *face2;
  face1 = (Face *) a; 
  face2 = (Face *) b;

  size_t compare1, compare2;

  if(face1->numvertices != face2->numvertices) {
    cerr << "Number of vertices mismatch in sortcompare!!";
    return 0;
  } 

  #ifdef DEBUG
  cout << "About to compare: " << face1->numvertices << " vertices in sortcompare";
  #endif

  for(int i = 0; i < face1->numvertices; i++) {

    compare1 = (size_t)(face1->vertices[i]);
    compare2 = (size_t)(face2->vertices[i]);

    #ifdef DEBUG
    cout << "\n" << compare1 << " vs " << compare2;
    #endif

    if(compare1 < compare2)
      return -1;
    if(compare1 > compare2)
      return 1;

  }

  return 0;

}

void Stack::Print() 
{

  if(totalblocks > 1) {
    cerr << "We can't dump big stacks yet!";
    exit(1);
  }

  for(int i = 0; i < totalsize; i++)
     ((Face *)((blocks[0])[i]))->PrintVertices(); 

  return;

}
