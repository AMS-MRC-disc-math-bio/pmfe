// Copyright (c) 2014 Andrew Gainer-Dewar.

#include <utility>
#include <gmpxx.h>
#include <cmath>
#include <iostream>

#include "iB4e.h"
#include "parametrizer_types.h"

std::ostream& operator<<(std::ostream& os, const ParameterVector& params) {
    os << "Multiloop penalty: " << params.multiloop_penalty.get_str(10) << std::endl
       << "Unpaired penalty: " << params.unpaired_penalty.get_str(10) << std::endl
       << "Branch penalty: " << params.branch_penalty.get_str(10) << std::endl
       << "Dummy scaling: " << params.dummy_scaling.get_str(10) << std::endl;
    return os;
};

ParameterVector::ParameterVector(QVector v) {
    assert(v.dimension() == 4);
    multiloop_penalty = mpq_class(v.cartesian(0).mpq());
    unpaired_penalty = mpq_class(v.cartesian(1).mpq());
    branch_penalty = mpq_class(v.cartesian(2).mpq());
    dummy_scaling = mpq_class(v.cartesian(3).mpq());
    this->canonicalize();
};

QVector ParameterVector::as_QVector() {
    this->canonicalize();
    mpq_t values [4];
    mpq_inits(values[0], values[1], values[2], values[3], NULL);

    mpq_set(values[0], multiloop_penalty.get_mpq_t());
    mpq_set(values[1], unpaired_penalty.get_mpq_t());
    mpq_set(values[2], branch_penalty.get_mpq_t());
    mpq_set(values[3], dummy_scaling.get_mpq_t());

    QVector result(4, values, values+4);
    mpq_clears(values[0], values[1], values[2], values[3], NULL);
    return result;
};

std::ostream& operator<<(std::ostream& os, const ScoreVector& score) {
    os << "Multiloops: " << score.multiloops.get_str(10) << std::endl
       << "Unpaired bases: " << score.unpaired.get_str(10) << std::endl
       << "Branches: " << score.branches.get_str(10) << std::endl
       << "w: " << score.w.get_str(10) << std::endl
       << "Parametrized energy: " << score.energy.get_str(10) << std::endl;
    return os;
};

QPoint ScoreVector::get_q4point() {
    this->canonicalize();
    mpq_t values [4];
    mpq_inits(values[0], values[1], values[2], values[3], NULL);

    mpq_set_z(values[0], multiloops.get_mpz_t());
    mpq_set_z(values[1], unpaired.get_mpz_t());
    mpq_set_z(values[2], branches.get_mpz_t());
    mpq_set(values[3], w.get_mpq_t());

    QPoint result(4, values, values+4);
    mpq_clears(values[0], values[1], values[2], values[3], NULL);
    return result;
};

std::ostream& operator<<(std::ostream& os, const energy_pair& energy) {
    os << "(" << energy.param.get_str(10) << ", " << energy.classical.get_str(10) << ")";
    return os;
};

bool operator< (const energy_pair &a, const energy_pair &b) {
    return (a.param < b.param);
};
bool operator<= (const energy_pair &a, const energy_pair &b) {
    return (a.param <= b.param);
};

bool operator== (const energy_pair &a, const energy_pair &b) {
    return (a.param == b.param);
};
bool operator!= (const energy_pair &a, const energy_pair &b) {
    return not (a == b);
};

bool operator> (const energy_pair &a, const energy_pair &b) {
    return (a.param > b.param);
};
bool operator>= (const energy_pair &a, const energy_pair &b) {
    return (a.param >= b.param);
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

mpq_class multiloop_default = mpq_class(34, 10);
mpq_class unpaired_default = mpq_class(0);
mpq_class branch_default = mpq_class(4, 10);
mpq_class dummy_default = mpq_class(1);
