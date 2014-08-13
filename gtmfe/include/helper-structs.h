#ifndef _HELPER_STRUCTS_H_
#define _HELPER_STRUCTS_H_

struct PolytopeVector {
  int multiloops, branches, unpaired;
  long double w, energy;
};

struct ParameterVector {
  long double a, b, c, d;
};

#endif
