#ifndef _HELPER_STRUCTS_H_
#define _HELPER_STRUCTS_H_

#include <vector>
#include <utility>
#include <gmpxx.h>

class ParameterVector {
 public:
  mpq_class multiloop_penalty, unpaired_penalty, branch_penalty, dummy_scaling;
 ParameterVector(mpq_class a = mpq_class(34,10), mpq_class b = mpq_class(0), mpq_class c = mpq_class(4,10), mpq_class d = mpq_class(1)) : multiloop_penalty(a), unpaired_penalty(b), branch_penalty(c), dummy_scaling(d) {};

  void set_from_pairs(std::pair<long, long> multiloop_pair, std::pair<long, long> unpaired_pair, std::pair<long, long> branch_pair, std::pair<long, long> dummy_pair);

  void set_from_words(const char* multiloop_word, const char* unpaired_word, const char* branch_word, const char* dummy_word);
  
  std::vector< std::pair<long, long> > get_pairs();
};

class PolytopeVector {
 public:
  int multiloops, branches, unpaired;
  mpq_class w, energy;
  //ParameterVector params;

  std::vector<std::pair<long, long> > get_pairs();
};

mpq_class get_mpq_from_word(const char* word);

#endif
