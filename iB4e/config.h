
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


#ifndef CONFIG_H
#define CONFIG_H


/*****************************************************************************/
  #define GMP_RATIONALS 100
  #define GMP_INTS 101
  #define LONGLONG_INTS 102

  /*******************************************************************
  **                                                                **
  **  Choose desired numerical data type, by defining NUMBER_TYPE   **
  **  to be GMP_RATIONALS, GMP_INTS, or LONGLONG_INTS               **
  **                                                                **
  ********************************************************************/
  #define NUMBER_TYPE LONGLONG_INTS



/*****************************************************************************/ 











  /****************************/

  //#define DEBUG

  //   #define DEBUG2
  //   #define DEBUGEXITS

  #define MAXSTACKBLOCKS 10000
  #define STACKBLOCKSIZE 500000

  #if NUMBER_TYPE == GMP_RATIONALS
    #define NUMBER mpq_class
    #define GMP
  #endif
  
  #if NUMBER_TYPE == GMP_INTS
    #define NUMBER mpz_class
    #define GMP
  #endif

  #if NUMBER_TYPE == LONGLONG_INTS
    #define NUMBER long long int
  #endif

  #ifdef GMP
    #include <gmpxx.h>
  #endif

  void printNumber(NUMBER a);


#endif
