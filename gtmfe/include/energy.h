#ifndef _ENERGY_TABLES_H_
#define _ENERGY_TABLES_H_

#include "data.h"
#include <gmpxx.h>

extern mpq_class *V; 
extern mpq_class *W; 
extern mpq_class *VBI; 
extern mpq_class *VM; 
extern mpq_class **WM; 
extern mpq_class **WMPrime; 
extern int *indx; 
extern mpq_class **PP; 

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
  mpq_class Ed3(int i, int j, int k);
  mpq_class Ed5(int i, int j, int k);
  mpq_class auPenalty(int i, int j);

  mpq_class eS(int i, int j);
  mpq_class eH(int i, int j);
  mpq_class eL(int i, int j, int ip, int jp);
  mpq_class eL1(int i, int j, int ip, int jp);
  mpq_class Estackm(int i, int j);
  mpq_class Estacke(int i, int j);

  void create_tables(int len);
  void init_tables(int len);
  void free_tables(int len);
#ifdef __cplusplus
}
#endif

#endif
