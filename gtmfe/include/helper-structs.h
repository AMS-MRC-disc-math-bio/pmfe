#ifndef _HELPER_STRUCTS_H_
#define _HELPER_STRUCTS_H_

#include <vector>
#include <utility>
#include <gmpxx.h>

class ParameterVector {
 public:
  mpq_class multiloop_penalty, unpaired_penalty, branch_penalty, dummy_scaling;
 ParameterVector(mpq_class a = mpq_class(34,10), mpq_class b = mpq_class(0), mpq_class c = mpq_class(4,10), mpq_class d = mpq_class(1)) : multiloop_penalty(a), unpaired_penalty(b), branch_penalty(c), dummy_scaling(d) {};

  void set_from_pairs(std::pair<long, long> multiloop_pair, std::pair<long, long> unpaired_pair, std::pair<long, long> branch_pair, std::pair<long, long> dummy_pair) {
    mpq_class multiloop_penalty(multiloop_pair.first, multiloop_pair.second);
    mpq_class unpaired_penalty(unpaired_pair.first, unpaired_pair.second);
    mpq_class branch_penalty(branch_pair.first, branch_pair.second);
    mpq_class dummy_scaling(dummy_pair.first, dummy_pair.second);
  };

  std::vector< std::pair<long, long> > get_pairs() {
    std::vector< std::pair<long, long> > result;
    result.push_back(std::pair<long, long> (multiloop_penalty.get_num().get_si(), multiloop_penalty.get_den().get_si()));
    result.push_back(std::pair<long, long> (unpaired_penalty.get_num().get_si(), unpaired_penalty.get_den().get_si()));
    result.push_back(std::pair<long, long> (branch_penalty.get_num().get_si(), branch_penalty.get_den().get_si()));
    result.push_back(std::pair<long, long> (dummy_scaling.get_num().get_si(), dummy_scaling.get_den().get_si()));
    return result;
  }
    
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
