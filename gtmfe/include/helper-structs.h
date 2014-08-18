#ifndef _HELPER_STRUCTS_H_
#define _HELPER_STRUCTS_H_

#include <gmpxx.h>

struct PolytopeVector {
  int multiloops, branches, unpaired;
  mpq_class w, energy;
};

struct ParameterVector {
  mpq_class a, b, c, d;
};

#endif
