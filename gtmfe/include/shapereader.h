#ifndef _SHAPEREADER_H_
#define _SHAPEREADER_H_

extern float* SHAPEarray;
extern int* SHAPEenergies;

void readSHAPEarray(const char* filename, int seqlength);
void free_shapeArray(int len);
void print_shapeArray(int len);
int calcShapeEnergy(float shapeNumber);
float shapeModel(float SHAPE_value);

#ifdef __cplusplus
extern "C"{
#endif

  int getShapeEnergy(int position);

#ifdef __cplusplus
}
#endif

#endif
