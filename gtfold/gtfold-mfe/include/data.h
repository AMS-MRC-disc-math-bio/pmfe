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

extern double poppen[5];
extern double maxpen;
extern double eparam[11];
extern double multConst[3]; /* for multiloop penalties. */
extern double dangle[4][4][4][2]; /* Contain dangling energy values */
extern double inter[31]; /* Contains size penalty for internal loops */
extern double bulge[31]; /* Contain the size penalty for bulges */
extern double hairpin[31]; /* Contains the size penalty for hairpin loops */
extern double stack[256]; /* Stacking energy used to calculate energy of stack loops */
extern double tstkh[256]; /* Terminal mismatch energy used in the calculations of hairpin loops */
extern double tstki[256]; /* Terminal mismatch energy used in the calculations of internal loops */
extern double tloop[maxtloop + 1][2];
extern int numoftloops;
extern double iloop22[5][5][5][5][5][5][5][5]; /* 2*2 internal looops */
extern double iloop21[5][5][5][5][5][5][5]; /* 2*1 internal loops */
extern double iloop11[5][5][5][5][5][5]; /* 1*1 internal loops */
extern double coax[6][6][6][6];
extern double tstackcoax[6][6][6][6];
extern double coaxstack[6][6][6][6];
extern double tstack[6][6][6][6];
extern double tstkm[6][6][6][6];
extern double auend;
extern double gubonus;
extern double cint; /* cint, cslope, c3 are used for poly C hairpin loops */
extern double cslope;
extern double c3;
extern double efn2a;
extern double efn2b;
extern double efn2c;
extern double triloop[maxtloop + 1][2];
extern int numoftriloops;
extern double init;
extern int gail; /* It is either 0 or 1. It is used for grosely asymmetric internal loops */
extern double prelog;

extern double tstackm[5][5][6][6];
extern double tstacke[5][5][6][6];
extern double tstacki23[5][5][5][5];


#define fourBaseIndex(a, b, c, d) (((a) << 6) + ((b) << 4) + ((c) << 2) + (d))

#endif

