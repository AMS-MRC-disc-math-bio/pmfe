
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


  /****************************/

  //#define DEBUG

  //   #define DEBUG2
  //   #define DEBUGEXITS

  #define MAXSTACKBLOCKS 10000
  #define STACKBLOCKSIZE 500000

  #include <gmpxx.h>

  void printNumber(mpq_class a);


#endif
