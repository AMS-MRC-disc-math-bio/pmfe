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

extern long double poppen[5];
extern long double maxpen;
extern long double eparam[11];
extern long double multConst[3]; /* for multiloop penalties. */
extern long double dangle[4][4][4][2]; /* Contain dangling energy values */
extern long double inter[31]; /* Contains size penalty for internal loops */
extern long double bulge[31]; /* Contain the size penalty for bulges */
extern long double hairpin[31]; /* Contains the size penalty for hairpin loops */
extern long double stack[256]; /* Stacking energy used to calculate energy of stack loops */
extern long double tstkh[256]; /* Terminal mismatch energy used in the calculations of hairpin loops */
extern long double tstki[256]; /* Terminal mismatch energy used in the calculations of internal loops */
extern long double tloop[maxtloop + 1][2];
extern int numoftloops;
extern long double iloop22[5][5][5][5][5][5][5][5]; /* 2*2 internal looops */
extern long double iloop21[5][5][5][5][5][5][5]; /* 2*1 internal loops */
extern long double iloop11[5][5][5][5][5][5]; /* 1*1 internal loops */
extern long double coax[6][6][6][6];
extern long double tstackcoax[6][6][6][6];
extern long double coaxstack[6][6][6][6];
extern long double tstack[6][6][6][6];
extern long double tstkm[6][6][6][6];
extern long double auend;
extern long double gubonus;
extern long double cint; /* cint, cslope, c3 are used for poly C hairpin loops */
extern long double cslope;
extern long double c3;
extern long double efn2a;
extern long double efn2b;
extern long double efn2c;
extern long double triloop[maxtloop + 1][2];
extern int numoftriloops;
extern long double init;
extern int gail; /* It is either 0 or 1. It is used for grosely asymmetric internal loops */
extern long double prelog;

extern long double tstackm[5][5][6][6];
extern long double tstacke[5][5][6][6];
extern long double tstacki23[5][5][5][5];


#define fourBaseIndex(a, b, c, d) (((a) << 6) + ((b) << 4) + ((c) << 2) + (d))

#endif

