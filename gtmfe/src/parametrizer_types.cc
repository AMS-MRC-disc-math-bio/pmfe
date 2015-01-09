// Copyright (c) 2014 Andrew Gainer-Dewar.

#include <utility>
#include <vector>
#include <gmpxx.h>
#include <cmath>
#include <iostream>

#include "parametrizer_types.h"

#include <boost/python/tuple.hpp>
#include <boost/python/extract.hpp>

namespace py = boost::python;

mpq_class mpq_from_pair(py::tuple pair) {
    mpz_class num(0), den(1);
    py::extract<int> numext(pair[0]), denext(pair[1]);

    if (numext.check() && denext.check())
    {
        num = numext();
        den = denext();
    }

    mpq_class result (num, den);

    return result;
};

mpz_class mpz_from_pair(py::tuple pair) {
    mpz_class num(0), den(1);
    py::extract<int> numext(pair[0]), denext(pair[1]);

    if (numext.check() && denext.check())
    {
        num = numext();
        den = denext();
    }

    // TODO: Deal with case that den != 1
    mpz_class result (num);
    return result;
};

py::tuple pair_from_mpq(mpq_class value) {
    // Warning: this is limited to long precision!
    py::tuple pair =
        py::make_tuple(value.get_num().get_si(), value.get_den().get_si());

    return pair;
};

py::tuple pair_from_mpz(mpz_class value) {
    // Warning: this is limited to long precision!
    py::tuple pair =
        py::make_tuple(value.get_si(), 1);

    return pair;
};

std::ostream& operator<<(std::ostream& os, const ParameterVector& params) {
    os << "Multiloop penalty: " << params.multiloop_penalty.get_str(10) << std::endl
       << "Unpaired penalty: " << params.unpaired_penalty.get_str(10) << std::endl
       << "Branch penalty: " << params.branch_penalty.get_str(10) << std::endl
       << "Dummy scaling: " << params.dummy_scaling.get_str(10) << std::endl;
    return os;
};

bool operator==(const ParameterVector& a, const ParameterVector& b) {
    return (
            a.multiloop_penalty == b.multiloop_penalty &&
            a.unpaired_penalty == b.unpaired_penalty &&
            a.branch_penalty == b.branch_penalty &&
            a.dummy_scaling == b.dummy_scaling
            );
};

bool operator!=(const ParameterVector& a, const ParameterVector& b) {
    return !(a == b);
};

py::tuple ParameterVector::as_pairs() {
    py::tuple pairs =
        py::make_tuple(
                       pair_from_mpq(multiloop_penalty),
                       pair_from_mpq(unpaired_penalty),
                       pair_from_mpq(branch_penalty),
                       pair_from_mpq(dummy_scaling)
                       );

    return pairs;
};

ParameterVector::ParameterVector(py::tuple p_multiloop_penalty, py::tuple p_unpaired_penalty, py::tuple p_branch_penalty, py::tuple p_dummy_scaling) {
    multiloop_penalty = mpq_from_pair(p_multiloop_penalty);
    branch_penalty = mpq_from_pair(p_unpaired_penalty);
    unpaired_penalty = mpq_from_pair(p_branch_penalty);
    dummy_scaling = mpq_from_pair(p_dummy_scaling);
    this->canonicalize();
};

std::ostream& operator<<(std::ostream& os, const ScoreVector& score) {
    os << "Multiloops: " << score.multiloops.get_str(10) << std::endl
       << "Unpaired bases: " << score.unpaired.get_str(10) << std::endl
       << "Branches: " << score.branches.get_str(10) << std::endl
       << "w: " << score.w.get_str(10) << std::endl
       << "Parametrized energy: " << score.energy.get_str(10) << std::endl;
    return os;
};

bool operator==(const ScoreVector& a, const ScoreVector& b) {
    return (
            a.multiloops == b.multiloops &&
            a.unpaired == b.unpaired &&
            a.branches == b.branches &&
            a.w == b.w &&
            a.energy == b.energy
            );
};

bool operator!=(const ScoreVector& a, const ScoreVector& b) {
    return !(a == b);
}

py::tuple ScoreVector::as_pairs() {
    py::tuple pairs =
        py::make_tuple(
                       pair_from_mpz(multiloops),
                       pair_from_mpq(unpaired),
                       pair_from_mpz(branches),
                       pair_from_mpq(w),
                       pair_from_mpq(energy)
                       );
    return pairs;
};

ScoreVector::ScoreVector(py::tuple p_multiloops, py::tuple p_unpaired, py::tuple p_branches, py::tuple p_w, py::tuple p_energy) {
    multiloops = mpz_from_pair(p_multiloops);
    unpaired = mpz_from_pair(p_unpaired);
    branches = mpz_from_pair(p_branches);
    w = mpq_from_pair(p_w);
    energy = mpq_from_pair(p_energy);
    this->canonicalize();
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
  result.canonicalize();
  return result;
};

mpq_class multiloop_default = mpq_class(17, 5);
mpq_class unpaired_default = mpq_class(0);
mpq_class branch_default = mpq_class(2, 5);
mpq_class dummy_default = mpq_class(1);
