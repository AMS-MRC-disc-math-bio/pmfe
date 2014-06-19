
// Copyright Peter Huggins 2006
// Released under the GNU GPL license

/*
    This file is a part of iB4e, which is free software;
    you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/



#ifndef LINALG_H
#define LINALG_H


#include "config.h"

#if NUMBER_TYPE != GMP_RATIONALS
  #error linalg.cc requires GMP rational numbers
#endif


int uppertriangular(NUMBER *matrix, int numrows, int numcols);
void printmatrix(NUMBER *matrix, int numrows, int numcols);
NUMBER diagonalproduct(NUMBER *matrix, int numrows); 


#endif
