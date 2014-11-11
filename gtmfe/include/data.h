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

extern energy_pair poppen[5];
extern energy_pair maxpen;
extern energy_pair eparam[11];
extern energy_pair multConst[3]; /* for multiloop penalties. */
extern energy_pair dangle[4][4][4][2]; /* Contain dangling energy values */
extern energy_pair inter[31]; /* Contains size penalty for internal loops */
extern energy_pair bulge[31]; /* Contain the size penalty for bulges */
extern energy_pair hairpin[31]; /* Contains the size penalty for hairpin loops */
extern energy_pair stack[256]; /* Stacking energy used to calculate energy of stack loops */
extern energy_pair tstkh[256]; /* Terminal mismatch energy used in the calculations of hairpin loops */
extern energy_pair tstki[256]; /* Terminal mismatch energy used in the calculations of internal loops */
extern int tloop[maxtloop + 1][2];
extern int numoftloops;
extern energy_pair iloop22[5][5][5][5][5][5][5][5]; /* 2*2 internal looops */
extern energy_pair iloop21[5][5][5][5][5][5][5]; /* 2*1 internal loops */
extern energy_pair iloop11[5][5][5][5][5][5]; /* 1*1 internal loops */
extern energy_pair coax[6][6][6][6];
extern energy_pair tstackcoax[6][6][6][6];
extern energy_pair coaxstack[6][6][6][6];
extern energy_pair tstack[6][6][6][6];
extern energy_pair tstkm[6][6][6][6];
extern energy_pair auend;
extern energy_pair gubonus;
extern energy_pair cint; /* cint, cslope, c3 are used for poly C hairpin loops */
extern energy_pair cslope;
extern energy_pair c3;
extern energy_pair efn2a;
extern energy_pair efn2b;
extern energy_pair efn2c;
extern energy_pair triloop[maxtloop + 1][2];
extern int numoftriloops;
extern energy_pair init;
extern int gail; /* It is either 0 or 1. It is used for grosely asymmetric internal loops */
extern energy_pair prelog;

extern energy_pair tstackm[5][5][6][6];
extern energy_pair tstacke[5][5][6][6];
extern energy_pair tstacki23[5][5][5][5];

extern energy_pair dummy_scaling;
extern energy_pair multiloop_penalty;
extern energy_pair unpaired_penalty;
extern energy_pair branch_penalty;


#define fourBaseIndex(a, b, c, d) (((a) << 6) + ((b) << 4) + ((c) << 2) + (d))

#endif
