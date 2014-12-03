// Copyright (c) 2014 Andrew Gainer-Dewar.

#ifndef PARAMETRIZER_TYPES_H
#define PARAMETRIZER_TYPES_H

#include <string>
#include <utility>
#include <gmpxx.h>
#include <iostream>

#include "iB4e.h"

extern mpq_class multiloop_default;
extern mpq_class unpaired_default;
extern mpq_class branch_default;
extern mpq_class dummy_default;

class ParameterVector {
 public:
 ParameterVector(mpq_class multiloop_penalty = multiloop_default, mpq_class unpaired_penalty = unpaired_default, mpq_class branch_penalty = branch_default, mpq_class dummy_scaling = dummy_default) : multiloop_penalty(multiloop_penalty), unpaired_penalty(unpaired_penalty), branch_penalty(branch_penalty), dummy_scaling(dummy_scaling) {
        this->canonicalize();
    };
    mpq_class multiloop_penalty, unpaired_penalty, branch_penalty, dummy_scaling;

    ParameterVector(QVector v);

    void canonicalize() {
        multiloop_penalty.canonicalize();
        unpaired_penalty.canonicalize();
        branch_penalty.canonicalize();
        dummy_scaling.canonicalize();
    }

    QVector as_QVector();

    friend std::ostream& operator<<(std::ostream& os, const ParameterVector& params);
};

class ScoreVector {
 public:
 ScoreVector(mpz_class multiloops = 0, mpz_class branches = 0, mpz_class unpaired = 0, mpq_class w = mpq_class(0), mpq_class energy = mpq_class(0)) : multiloops(multiloops), branches(branches), unpaired(unpaired), w(w), energy(energy) {
        this->canonicalize();
    };
    mpz_class multiloops, branches, unpaired;
    mpq_class w, energy;

    QPoint get_q4point();

    void canonicalize(){
        w.canonicalize();
        energy.canonicalize();
    }

    friend std::ostream& operator<<(std::ostream& os, const ScoreVector& score);
};

mpq_class get_mpq_from_word(std::string word);

#endif
