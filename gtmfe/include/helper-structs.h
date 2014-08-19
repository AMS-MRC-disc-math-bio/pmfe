#ifndef _HELPER_STRUCTS_H_
#define _HELPER_STRUCTS_H_

#include <vector>
#include <utility>
#include <gmpxx.h>

class ParameterVector {
 public:
  mpq_class multiloop_penalty, unpaired_penalty, branch_penalty, dummy_scaling;
 ParameterVector(mpq_class a = mpq_class(34,10), mpq_class b = mpq_class(0), mpq_class c = mpq_class(4,10), mpq_class d = mpq_class(1)) : multiloop_penalty(a), unpaired_penalty(b), branch_penalty(c), dummy_scaling(d) {};
};

class PolytopeVector {
 public:
  int multiloops, branches, unpaired;
  mpq_class w, energy;
  ParameterVector params;

  std::vector< std::pair<long, long> > get_pairs() {
    std::vector< std::pair<long, long> > result;
    result.push_back(std::pair<long, long> (multiloops, 1));
    result.push_back(std::pair<long, long> (branches, 1));
    result.push_back(std::pair<long, long> (unpaired, 1));
    result.push_back(std::pair<long, long> (w.get_num().get_si(), w.get_den().get_si()));
    result.push_back(std::pair<long, long> (energy.get_num().get_si(), energy.get_den().get_si()));
    return result;
  }
};

#endif
