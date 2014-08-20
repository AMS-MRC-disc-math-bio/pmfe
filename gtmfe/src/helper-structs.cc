#include "helper-structs.h"
#include "constants.h"
#include <vector>
#include <utility>
#include <gmpxx.h>
#include <cmath>

void ParameterVector::set_from_pairs(std::pair<long, long> multiloop_pair, std::pair<long, long> unpaired_pair, std::pair<long, long> branch_pair, std::pair<long, long> dummy_pair) {
  multiloop_penalty = mpq_class (multiloop_pair.first, multiloop_pair.second);
  multiloop_penalty.canonicalize();
    
  unpaired_penalty = mpq_class (unpaired_pair.first, unpaired_pair.second);
  unpaired_penalty.canonicalize();
    
  branch_penalty = mpq_class (branch_pair.first, branch_pair.second);
  branch_penalty.canonicalize();
    
  dummy_scaling = mpq_class (dummy_pair.first, dummy_pair.second);
  dummy_scaling.canonicalize();
};

void ParameterVector::set_from_words(const char* multiloop_word, const char* unpaired_word, const char* branch_word, const char* dummy_word) {
  multiloop_penalty = get_mpq_from_word(multiloop_word);
  unpaired_penalty = get_mpq_from_word(unpaired_word);
  branch_penalty = get_mpq_from_word(branch_word);
  dummy_scaling = get_mpq_from_word(dummy_word);
};
  

std::vector< std::pair<long, long> > ParameterVector::get_pairs() {
  std::vector< std::pair<long, long> > result;
  result.push_back(std::pair<long, long> (multiloop_penalty.get_num().get_si(), multiloop_penalty.get_den().get_si()));
  result.push_back(std::pair<long, long> (unpaired_penalty.get_num().get_si(), unpaired_penalty.get_den().get_si()));
  result.push_back(std::pair<long, long> (branch_penalty.get_num().get_si(), branch_penalty.get_den().get_si()));
  result.push_back(std::pair<long, long> (dummy_scaling.get_num().get_si(), dummy_scaling.get_den().get_si()));
  return result;
};

std::vector<std::pair<long, long> > PolytopeVector::get_pairs() {
  std::vector< std::pair<long, long> > result;
  result.push_back(std::pair<long, long> (multiloops, 1));
  result.push_back(std::pair<long, long> (unpaired, 1));
  result.push_back(std::pair<long, long> (branches, 1));
  result.push_back(std::pair<long, long> (w.get_num().get_si(), w.get_den().get_si()));
  result.push_back(std::pair<long, long> (energy.get_num().get_si(), energy.get_den().get_si()));
  return result;
};

mpq_class get_mpq_from_word(const char* word) {
  mpq_class result;
  std::string decimalword (word);
  std::size_t decimalpoint = decimalword.find('.');
  bool negative = (decimalword.find('-') != std::string::npos);
  if (decimalpoint != std::string::npos) {
    std::string intpart = decimalword.substr(0, decimalpoint);
    std::string fracpart = decimalword.substr(decimalpoint+1);

    if (intpart == "") intpart = "0";
    mpz_class theint (intpart, 10);
    theint = abs(theint);

    // Carve out the fractional part. Surprisingly fiddly!
    int fracdenom = pow(10, fracpart.length());
    mpz_class fracval (fracpart, 10);
    mpq_class thefrac (fracval, fracdenom);
    thefrac.canonicalize();

    result = theint + thefrac;
    if (negative) result *= -1;
  } else {
    mpq_class thevalue (decimalword);
    result = thevalue;
  }
  return result;
};
