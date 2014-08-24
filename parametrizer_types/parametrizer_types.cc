#include "parametrizer_types.h"
#include <vector>
#include <utility>
#include <gmpxx.h>
#include <cmath>

void ParameterVector::set_from_pairs(std::pair<long, long> multiloop_param_pair, std::pair<long, long> unpaired_param_pair, std::pair<long, long> branch_param_pair, std::pair<long, long> dummy_param_pair) {
  multiloop_penalty = mpq_class (multiloop_param_pair.first, multiloop_param_pair.second);
  multiloop_penalty.canonicalize();
    
  unpaired_penalty = mpq_class (unpaired_param_pair.first, unpaired_param_pair.second);
  unpaired_penalty.canonicalize();
    
  branch_penalty = mpq_class (branch_param_pair.first, branch_param_pair.second);
  branch_penalty.canonicalize();
    
  dummy_scaling = mpq_class (dummy_param_pair.first, dummy_param_pair.second);
  dummy_scaling.canonicalize();
};

void ParameterVector::set_from_words(std::string multiloop_param_word, std::string unpaired_param_word, std::string branch_param_word, std::string dummy_param_word) {
  multiloop_penalty = get_mpq_from_word(multiloop_param_word);
  unpaired_penalty = get_mpq_from_word(unpaired_param_word);
  branch_penalty = get_mpq_from_word(branch_param_word);
  dummy_scaling = get_mpq_from_word(dummy_param_word);
};

std::vector< std::pair<long, long> > ParameterVector::get_pairs() {
  std::vector< std::pair<long, long> > result;
  result.push_back(std::pair<long, long> (multiloop_penalty.get_num().get_si(), multiloop_penalty.get_den().get_si()));
  result.push_back(std::pair<long, long> (unpaired_penalty.get_num().get_si(), unpaired_penalty.get_den().get_si()));
  result.push_back(std::pair<long, long> (branch_penalty.get_num().get_si(), branch_penalty.get_den().get_si()));
  result.push_back(std::pair<long, long> (dummy_scaling.get_num().get_si(), dummy_scaling.get_den().get_si()));
  return result;
};

std::vector< std::string > ParameterVector::get_words() {
  std::vector< std::string > result;
  result.push_back(multiloop_penalty.get_str(10));
  result.push_back(unpaired_penalty.get_str(10));
  result.push_back(branch_penalty.get_str(10));
  result.push_back(dummy_scaling.get_str(10));
  return result;
};

void ScoreVector::set_from_pairs(long multiloop_score, long unpaired_score, long branch_score, std::pair<long, long> w_score_pair, std::pair<long, long> energy_pair) {
  multiloops = multiloop_score;
  unpaired = unpaired_score;
  branches = branch_score;
  
  w = mpq_class (w_score_pair.first, w_score_pair.second);
  w.canonicalize();
  
  energy = mpq_class (energy_pair.first, energy_pair.second);
  energy.canonicalize();
};

void ScoreVector::set_from_words(std::string multiloop_score_word, std::string unpaired_score_word, std::string branch_score_word, std::string w_score_word, std::string energy_word) {
  multiloops.set_str(multiloop_score_word, 10);
  unpaired.set_str(unpaired_score_word, 10);
  branches.set_str(branch_score_word, 10);
  w = get_mpq_from_word(w_score_word);
  energy = get_mpq_from_word(energy_word);
};

std::vector< std::pair<long, long> > ScoreVector::get_pairs() {
  std::vector< std::pair<long, long> > result;
  result.push_back(std::pair<long, long> (multiloops.get_si(), 1));
  result.push_back(std::pair<long, long> (unpaired.get_si(), 1));
  result.push_back(std::pair<long, long> (branches.get_si(), 1));
  result.push_back(std::pair<long, long> (w.get_num().get_si(), w.get_den().get_si()));
  result.push_back(std::pair<long, long> (energy.get_num().get_si(), energy.get_den().get_si()));
  return result;
};

std::vector< std::string > ScoreVector::get_words() {
  std::vector< std::string > result;
  result.push_back(multiloops.get_str(10));
  result.push_back(unpaired.get_str(10));
  result.push_back(branches.get_str(10));
  result.push_back(w.get_str(10));
  result.push_back(energy.get_str(10));
  return result;
};

mpq_class get_mpq_from_word(std::string word) {
  mpq_class result;
  std::size_t decimalpoint = word.find('.');
  bool negative = (word.find('-') != std::string::npos);
  if (decimalpoint != std::string::npos) {
    std::string intpart = word.substr(0, decimalpoint);
    std::string fracpart = word.substr(decimalpoint+1);

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
    mpq_class thevalue (word);
    result = thevalue;
  }
  return result;
};
