#ifndef _HELPER_STRUCTS_H_
#define _HELPER_STRUCTS_H_

struct PolytopeVector {
  int multiloops, branches, unpaired;
  double w, energy;
};

struct ParameterVector {
  double a, b, c, d;
};

#endif
