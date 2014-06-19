// Copyright 2006 Peter Huggins
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



#include "config.h"

#if NUMBER_TYPE == GMP_RATIONALS

#include "linalg.h"
#include <gmpxx.h>
#include <iostream>

// computes uppertriangular form of a matrix via unimodular transformation
// returns rank

using namespace std;


void printmatrix(NUMBER *matrix, int m, int n)
{
  for(int i = 0; i < m; i++) { 
    for(int j = 0; j < n; j++) { 
      printNumber(matrix[i*n + j]);
      cout << " ";
    }
    cout << "\n";
  }
  return;
}


int uppertriangular(NUMBER *matrix, int m, int n)
{

  int rank = 0;
  NUMBER swap;

  //int minmn;
  //if(m < n)
  //  minmn = m;
  //else
  //  minmn = n;

  int i = 0;
  for(int j = 0; j < n && i < m; j++) {
    // find nonzero entry in jth column (not above ith row)
    int nextrow;
    for(nextrow = i; nextrow < m; nextrow++)
      if(matrix[nextrow*n + j] != 0)
        break;

    if(nextrow < m) {
      rank++;
      if(nextrow != i) {
        // swap/negate row w/ nonzero entry with the ith row
        for(int j2 = 0; j2 < n; j2++) {
          swap = matrix[nextrow*n + j2];
          matrix[nextrow*n + j2] = matrix[i*n + j2];
          matrix[i*n + j2] = -1 * swap;
        }
      }

      // subtract appropriate multiples of ith row from all lower rows
      NUMBER multiple;
      for(int i2 = i+1; i2 < m; i2++) {
        multiple = matrix[i2*n + j] / matrix[i*n + j];
        for(int j2=0; j2 < n; j2++)
          matrix[i2*n + j2] -= multiple * matrix[i*n + j2];
      }
      i++;
    }

  }

  return rank;

}
      
NUMBER diagonalproduct(NUMBER *matrix, int n)
{
  NUMBER ans;
  ans = 1;
  for(int i = 0; i < n; i++)
    ans *= matrix[i*n + i];
  return ans;
}

#endif
