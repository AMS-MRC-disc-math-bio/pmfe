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

extern float poppen[5];
extern float maxpen;
extern float eparam[11];
extern float multConst[3]; /* for multiloop penalties. */
extern float dangle[4][4][4][2]; /* Contain dangling energy values */
extern float inter[31]; /* Contains size penalty for internal loops */
extern float bulge[31]; /* Contain the size penalty for bulges */
extern float hairpin[31]; /* Contains the size penalty for hairpin loops */
extern float stack[256]; /* Stacking energy used to calculate energy of stack loops */
extern float tstkh[256]; /* Terminal mismatch energy used in the calculations of hairpin loops */
extern float tstki[256]; /* Terminal mismatch energy used in the calculations of internal loops */
extern float tloop[maxtloop + 1][2];
extern int numoftloops;
extern float iloop22[5][5][5][5][5][5][5][5]; /* 2*2 internal looops */
extern float iloop21[5][5][5][5][5][5][5]; /* 2*1 internal loops */
extern float iloop11[5][5][5][5][5][5]; /* 1*1 internal loops */
extern float coax[6][6][6][6];
extern float tstackcoax[6][6][6][6];
extern float coaxstack[6][6][6][6];
extern float tstack[6][6][6][6];
extern float tstkm[6][6][6][6];
extern float auend;
extern float gubonus;
extern float cint; /* cint, cslope, c3 are used for poly C hairpin loops */
extern float cslope;
extern float c3;
extern float efn2a;
extern float efn2b;
extern float efn2c;
extern float triloop[maxtloop + 1][2];
extern int numoftriloops;
extern float init;
extern int gail; /* It is either 0 or 1. It is used for grosely asymmetric internal loops */
extern float prelog;

extern float tstackm[5][5][6][6];
extern float tstacke[5][5][6][6];
extern float tstacki23[5][5][5][5];


#define fourBaseIndex(a, b, c, d) (((a) << 6) + ((b) << 4) + ((c) << 2) + (d))

#endif

