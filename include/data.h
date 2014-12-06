/*
  GTfold: compute minimum free energy of RNA secondary structure
  Copyright (C) 2008  David A. Bader
  http://www.cc.gatech.edu/~bader

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* Amrita: some of the variables defined in this file are not uesed. */

#ifndef _DATA_H
#define _DATA_H

#include "constants.h"
#include "parametrizer_types.h"
#include <gmpxx.h>

extern mpq_class poppen[5];
extern mpq_class maxpen;
extern mpq_class eparam[11];
extern mpq_class multConst[3]; /* for multiloop penalties. */
extern mpq_class dangle[4][4][4][2]; /* Contain dangling energy values */
extern mpq_class inter[31]; /* Contains size penalty for internal loops */
extern mpq_class bulge[31]; /* Contain the size penalty for bulges */
extern mpq_class hairpin[31]; /* Contains the size penalty for hairpin loops */
extern mpq_class stack[256]; /* Stacking energy used to calculate energy of stack loops */
extern mpq_class tstkh[256]; /* Terminal mismatch energy used in the calculations of hairpin loops */
extern mpq_class tstki[256]; /* Terminal mismatch energy used in the calculations of internal loops */
extern mpq_class tloop[maxtloop + 1][2];
extern int numoftloops;
extern mpq_class iloop22[5][5][5][5][5][5][5][5]; /* 2*2 internal looops */
extern mpq_class iloop21[5][5][5][5][5][5][5]; /* 2*1 internal loops */
extern mpq_class iloop11[5][5][5][5][5][5]; /* 1*1 internal loops */
extern mpq_class coax[6][6][6][6];
extern mpq_class tstackcoax[6][6][6][6];
extern mpq_class coaxstack[6][6][6][6];
extern mpq_class tstack[6][6][6][6];
extern mpq_class tstkm[6][6][6][6];
extern mpq_class auend;
extern mpq_class gubonus;
extern mpq_class cint; /* cint, cslope, c3 are used for poly C hairpin loops */
extern mpq_class cslope;
extern mpq_class c3;
extern mpq_class efn2a;
extern mpq_class efn2b;
extern mpq_class efn2c;
extern mpq_class triloop[maxtloop + 1][2];
extern int numoftriloops;
extern mpq_class init;
extern int gail; /* It is either 0 or 1. It is used for grosely asymmetric internal loops */
extern mpq_class prelog;

extern mpq_class tstackm[5][5][6][6];
extern mpq_class tstacke[5][5][6][6];
extern mpq_class tstacki23[5][5][5][5];

extern mpq_class dummy_scaling;
extern mpq_class multiloop_penalty;
extern mpq_class unpaired_penalty;
extern mpq_class branch_penalty;


#define fourBaseIndex(a, b, c, d) (((a) << 6) + ((b) << 4) + ((c) << 2) + (d))

#endif
