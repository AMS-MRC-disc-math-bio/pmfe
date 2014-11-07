#ifndef _PARAMETRIZER_TYPES_H_
#define _PARAMETRIZER_TYPES_H_

#include <vector>
#include <string>
#include <utility>
#include <gmpxx.h>

class ParameterVector {
 public:
 ParameterVector(mpq_class multiloop_penalty = mpq_class(34,10), mpq_class unpaired_penalty = mpq_class(0), mpq_class branch_penalty = mpq_class(4,10), mpq_class dummy_scaling = mpq_class(1)) : multiloop_penalty(multiloop_penalty), unpaired_penalty(unpaired_penalty), branch_penalty(branch_penalty), dummy_scaling(dummy_scaling) {
        multiloop_penalty.canonicalize();
        unpaired_penalty.canonicalize();
        branch_penalty.canonicalize();
        dummy_scaling.canonicalize();
    };
    mpq_class multiloop_penalty, unpaired_penalty, branch_penalty, dummy_scaling;

    void set_from_pairs(std::pair<long, long> multiloop_param_pair, std::pair<long, long> unpaired_param_pair, std::pair<long, long> branch_param_pair, std::pair<long, long> dummy_param_pair);
    void set_from_words(std::string multiloop_param_word, std::string unpaired_param_word, std::string branch_param_word, std::string dummy_param_word);

    std::vector< std::pair<long, long> > get_pairs();
    std::vector< std::string > get_words();
};

class ScoreVector {
 public:
 ScoreVector(mpz_class multiloops = 0, mpz_class branches = 0, mpz_class unpaired = 0, mpq_class w = mpq_class(0), mpq_class energy = mpq_class(0)) : multiloops(multiloops), branches(branches), unpaired(unpaired), w(w), energy(energy) {
        w.canonicalize();
        energy.canonicalize();
    };
    mpz_class multiloops, branches, unpaired;
    mpq_class w, energy;

    void set_from_pairs(long multiloop_score, long unpaired_score, long branch_score, std::pair<long, long> w_score_pair,std::pair<long, long> energy_pair);
    void set_from_words(std::string multiloop_score_word, std::string unpaired_score_word, std::string branch_score_word, std::string w_score_word, std::string energy_word);

    std::vector< std::pair<long, long> > get_pairs();
    std::vector< std::string > get_words();
};

mpq_class get_mpq_from_word(std::string word);

#endif
