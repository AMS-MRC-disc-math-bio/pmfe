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

/* Modified by Prashant Gaurav <pgaurav3@gatech.edu>, 09/07/2010 */
/* Fixed the incorrect reporting of multiloop energy */

#ifndef _TRACEBACK_H
#define _TRACEBACK_H

#include "helper-structs.h"

PolytopeVector trace(int len);

void traceW(int i);
float traceV(int i, int j);
float traceVM(int i, int j);
float traceVBI(int i, int j);
float traceWM(int i, int j);
float traceWMPrime(int i, int j);

#endif
