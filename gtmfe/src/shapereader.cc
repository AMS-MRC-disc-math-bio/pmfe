#include<iostream>
#include<fstream>
#include<string.h>
#include<stdio.h>
#include <stdlib.h>
#include <math.h>

#include "shapereader.h"
#include "global.h"

using namespace std;

long double* SHAPEarray;
int* SHAPEenergies;

//static bool SHAPE_ENABLED = false;

//extern int SHAPE_ENABLED;

/*
  void enable_shape(bool b) {
  SHAPE_ENABLED = b;
  }*/

void enable_shape(int b) {
  SHAPE_ENABLED = b;
}



void free_shapeArray(int len){
  free(SHAPEarray);
}

void print_shapeArray(int len){
  for(int i=0; i<len; i++){
    printf("SHAPE value for position %d: %f, energy contribution: %f kcal/mol \n", i, SHAPEarray[i], (long double)SHAPEenergies[i]/100);
  }

}

void readSHAPEarray(const char* filename, int seqlength){

  ifstream infile(filename);
  SHAPEarray = (long double*)malloc(sizeof(long double)*(seqlength+1));	
  SHAPEenergies = (int*)malloc(sizeof(int)*(seqlength+1));

  for(int i = 0; i<seqlength; i++){
    SHAPEarray[i] = -999;	
    SHAPEenergies[i] = 0;
  }

  string line;
  int position;
  long double  SHAPEnumber;
	
  while(getline(infile,line)>0){
    if(sscanf(line.c_str(), "%d %lf", &position, &SHAPEnumber)==2){
      if(position <= seqlength){
        SHAPEarray[position] = SHAPEnumber;
        SHAPEenergies[position] = calcShapeEnergy(SHAPEnumber);
      }
      else{
        printf("Invalid SHAPE position indicator (ignoring line): %s\n", line.c_str()); 
      }	
    }else{
      printf("Invalid line (ignoring): %s\n", line.c_str());
    }
  }

}

int getShapeEnergy(int position){
  if(SHAPE_ENABLED){
    return SHAPEenergies[position];
  }	
  else{
    return 0;
  }
}

int calcShapeEnergy(long double shapeNumber){
  //ZS: This function returns the free energy contribution as an integer. 
  if(shapeNumber<(long double)0){
    return 0;
  }
  else{
    long double energy = shapeModel(shapeNumber);
    return (int)floor(100.0*energy+ .5);
  }
}


long double shapeModel(long double SHAPE_value){
  //ZS: This function calculates the free energy contribution due to SHAPE. 
  long double m = 2.6;
  long double b = -0.8;
  return m*log(SHAPE_value+1)+b;
}




