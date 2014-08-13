#ifndef _ENERGY_TABLES_H_
#define _ENERGY_TABLES_H_

#include "data.h"

extern long double *V; 
extern long double *W; 
extern long double *VBI; 
extern long double *VM; 
extern long double **WM; 
extern long double **WMPrime; 
extern int *indx; 
extern long double **PP; 

#define V(i,j) V[indx[j]+i]
#define VM(i,j) VM[indx[j]+i]
#define WM(i,j) WM[i][j]
#define WMPrime(i,j) WMPrime[i][j]
#define WMU(i,j) WM[i][j]
#define WML(i,j) WM[j][i]
#define VBI(i,j) VBI[indx[j]+i]
//#define RT ((0.00198721 * 310.15) * 100.00)
extern const float RT;
extern const float RT_;

#define auPen(i, j) ((( (i)==BASE_U || (j)==BASE_U ) && ( (i)==BASE_A || (i)==BASE_G || (j)==BASE_A || (j)==BASE_G )) ? auend : 0)

#ifdef __cplusplus
extern "C" {
#endif
  long double Ed3(int i, int j, int k);
  long double Ed5(int i, int j, int k);
  long double auPenalty(int i, int j);

#define Ec multConst[1]
#define Eb multConst[2]
#define Ea multConst[0] 

  long double eS(int i, int j);
  long double eH(int i, int j);
  long double eL(int i, int j, int ip, int jp);
  long double eL1(int i, int j, int ip, int jp);
  long double Estackm(int i, int j);
  long double Estacke(int i, int j);

  void create_tables(int len);
  void init_tables(int len);
  void free_tables(int len);
#ifdef __cplusplus
}
#endif

#endif
