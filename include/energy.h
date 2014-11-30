#ifndef _ENERGY_TABLES_H_
#define _ENERGY_TABLES_H_

#include "data.h"
#include "loader.h"
#include "parametrizer_types.h"
#include <gmpxx.h>
#include <vector>

extern std::vector<energy_pair> V;
extern std::vector<energy_pair> W;
extern std::vector<energy_pair> VBI;
extern std::vector<energy_pair> VM;
extern std::vector< std::vector<energy_pair> > WM;
extern std::vector< std::vector<energy_pair> > WMPrime;

extern std::vector< std::vector<mpq_class> > PP;
extern int *indx;

#define V_f(i,j) V[indx[j]+i]
#define VM_f(i,j) VM[indx[j]+i]
#define WM_f(i,j) WM[i][j]
#define WMPrime_f(i,j) WMPrime[i][j]
#define WMU_f(i,j) WM[i][j]
#define WML_f(i,j) WM[j][i]
#define VBI_f(i,j) VBI[indx[j]+i]

//#define RT ((0.00198721 * 310.15) * 100.00)
extern const float RT;
extern const float RT_;

energy_pair Ed3(int i, int j, int k);
energy_pair Ed5(int i, int j, int k);
energy_pair auPenalty(int i, int j);

energy_pair eS(int i, int j);
energy_pair eH(int i, int j);
energy_pair eL(int i, int j, int ip, int jp);
energy_pair eL1(int i, int j, int ip, int jp);
energy_pair Estackm(int i, int j);
energy_pair Estacke(int i, int j);

void create_tables(int len);
void init_tables(int len);
void free_tables(int len);

#endif
