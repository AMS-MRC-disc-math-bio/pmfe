// Copyright (c) 2014 Andrew Gainer-Dewar.

#ifndef PARAMETRIZER_TYPES_H
#define PARAMETRIZER_TYPES_H

#include <vector>
#include <string>
#include <utility>
#include <gmpxx.h>
#include <iostream>

#include "iB4e.h"

extern mpq_class multiloop_default;
extern mpq_class unpaired_default;
extern mpq_class branch_default;

class ParameterVector {
 public:
 ParameterVector(mpq_class multiloop_penalty = multiloop_default, mpq_class unpaired_penalty = unpaired_default, mpq_class branch_penalty = branch_default, mpq_class dummy_scaling = mpq_class(1)) : multiloop_penalty(multiloop_penalty), unpaired_penalty(unpaired_penalty), branch_penalty(branch_penalty), dummy_scaling(dummy_scaling) {
        multiloop_penalty.canonicalize();
        unpaired_penalty.canonicalize();
        branch_penalty.canonicalize();
        dummy_scaling.canonicalize();
    };
    mpq_class multiloop_penalty, unpaired_penalty, branch_penalty, dummy_scaling;

    ParameterVector(QVector v);

    void set_from_pairs(std::pair<long, long> multiloop_param_pair, std::pair<long, long> unpaired_param_pair, std::pair<long, long> branch_param_pair, std::pair<long, long> dummy_param_pair);
    void set_from_words(std::string multiloop_param_word, std::string unpaired_param_word, std::string branch_param_word, std::string dummy_param_word);

    std::vector< std::pair<long, long> > get_pairs();
    std::vector< std::string > get_words();

    friend std::ostream& operator<<(std::ostream& os, const ParameterVector& params);
};

class ScoreVector {
 public:
 ScoreVector(mpz_class multiloops = 0, mpz_class branches = 0, mpz_class unpaired = 0, mpq_class w = mpq_class(0), mpq_class energy = mpq_class(0)) : multiloops(multiloops), branches(branches), unpaired(unpaired), w(w), energy(energy) {
        w.canonicalize();
        energy.canonicalize();
    };
    mpz_class multiloops, branches, unpaired;
    mpq_class w, energy;

    QPoint get_q4point();

    void set_from_pairs(long multiloop_score, long unpaired_score, long branch_score, std::pair<long, long> w_score_pair,std::pair<long, long> energy_pair);
    void set_from_words(std::string multiloop_score_word, std::string unpaired_score_word, std::string branch_score_word, std::string w_score_word, std::string energy_word);

    std::vector< std::pair<long, long> > get_pairs();
    std::vector< std::string > get_words();
    friend std::ostream& operator<<(std::ostream& os, const ScoreVector& score);
};

class energy_pair {
 public:
 energy_pair(mpq_class param = mpq_class(0), mpq_class classical = mpq_class(0)): param(param), classical(classical) {};

    mpq_class param, classical;

    energy_pair& operator= (const energy_pair b) {
        this->param = b.param;
        this->classical = b.classical;
        return *this;
    }

    energy_pair& operator+= (const energy_pair &b) {
        this->param += b.param;
        this->classical += b.classical;
        return *this;
    }
    friend energy_pair operator+ (energy_pair a, const energy_pair &b){
        return a += b;
    }

    energy_pair& operator-= (const energy_pair &b){
        this->param -= b.param;
        this->classical -= b.classical;
        return *this;
    }

    friend energy_pair operator- (energy_pair a, const energy_pair &b) {
        return a -= b;
    }

    energy_pair& operator*= (const energy_pair &b){
        this->param *= b.param;
        this->classical *= b.classical;
        return *this;
    }

    friend energy_pair operator* (energy_pair a, const energy_pair &b) {
        return a *= b;
    }

    energy_pair& operator*= (const mpq_class &b){
        return *this *= energy_pair(b, b);
    }

    friend energy_pair operator* (energy_pair a, const mpq_class &b){
        return a *= b;
    }

    friend energy_pair operator* (const mpq_class &a, energy_pair b) {
        return b *= a;
    }

    friend bool operator< (const energy_pair &a, const energy_pair &b) { return a.param < b.param; };
    friend bool operator<= (const energy_pair &a, const energy_pair &b) { return a.param <= b.param; };

    friend bool operator== (const energy_pair &a, const energy_pair &b) { return a.param == b.param; };
    friend bool operator!= (const energy_pair &a, const energy_pair &b) { return a.param != b.param; };

    friend bool operator> (const energy_pair &a, const energy_pair &b) { return a.param > b.param; };
    friend bool operator>= (const energy_pair &a, const energy_pair &b) { return a.param >= b.param; };
};


mpq_class get_mpq_from_word(std::string word);

#endif
