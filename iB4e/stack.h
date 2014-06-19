
#ifndef STACK_H
#define STACK_H

#include "config.h"


int sortcompare(void *a, void *b);


class Stack {

 public:
  
  void **blocks[MAXSTACKBLOCKS];
  int currentblock;
  int currentpos;
  int totalblocks;
  int totalsize;
  bool empty;

  Stack() { blocks[0] = new void *[STACKBLOCKSIZE]; 
            currentblock = 0; 
            currentpos = 0; 
            totalblocks = 0;
            empty = true;
            totalsize=0;
          };

  ~Stack() { for(int i=0; i<totalblocks; i++) delete blocks[i]; };

   void Push(void *);
   void * Pop();
   void Sort();
   void Print();


};


#endif
