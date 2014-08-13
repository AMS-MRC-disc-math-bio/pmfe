#ifndef _HELPER_STRUCTS_H_
#define _HELPER_STRUCTS_H_

struct PolytopeVector {
  int multiloops, branches, unpaired;
  float w, energy;
};

struct ParameterVector {
  float a, b, c, d;
};

#endif
