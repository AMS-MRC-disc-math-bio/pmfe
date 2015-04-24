// Copyright (c) 2014 Andrew Gainer-Dewar.

#ifndef PARAMETRIZER_TYPES_H
#define PARAMETRIZER_TYPES_H

#include <string>
#include <utility>
#include <gmpxx.h>
#include <iostream>

#include <deque>

#include <boost/python/tuple.hpp>
#include <boost/filesystem.hpp>

#include "interval_tree.h"

namespace pmfe {
    extern const mpq_class multiloop_default;
    extern const mpq_class unpaired_default;
    extern const mpq_class branch_default;
    extern const mpq_class dummy_default;

    namespace py = boost::python;
    namespace fs = boost::filesystem;

    class ParameterVector {
    public:
    ParameterVector(mpq_class multiloop_penalty = multiloop_default, mpq_class unpaired_penalty = unpaired_default, mpq_class branch_penalty = branch_default, mpq_class dummy_scaling = dummy_default) : multiloop_penalty(multiloop_penalty), unpaired_penalty(unpaired_penalty), branch_penalty(branch_penalty), dummy_scaling(dummy_scaling) {
            this->canonicalize();
        };
        mpq_class multiloop_penalty, unpaired_penalty, branch_penalty, dummy_scaling;

        ParameterVector(py::tuple p_multiloop_penalty, py::tuple p_unpaired_penalty, py::tuple p_branch_penalty, py::tuple p_dummy_scaling);

        void canonicalize() {
            multiloop_penalty.canonicalize();
            unpaired_penalty.canonicalize();
            branch_penalty.canonicalize();
            dummy_scaling.canonicalize();
        }

        boost::python::tuple as_pairs();

        std::string print_as_list();

        friend std::ostream& operator<<(std::ostream& os, const ParameterVector& params);
        friend bool operator==(const ParameterVector& a, const ParameterVector& b);
        friend bool operator!=(const ParameterVector& a, const ParameterVector& b);
    };

    class ScoreVector {
    public:
    ScoreVector(mpz_class multiloops = 0, mpz_class unpaired = 0, mpz_class branches = 0, mpq_class w = mpq_class(0), mpq_class energy = mpq_class(0)) : multiloops(multiloops), branches(branches), unpaired(unpaired), w(w), energy(energy) {
            this->canonicalize();
        };
        mpz_class multiloops, branches, unpaired;
        mpq_class w, energy;

        ScoreVector(py::tuple p_multiloops, py::tuple p_unpaired, py::tuple p_branches, py::tuple p_w, py::tuple p_energy);

        boost::python::tuple as_pairs();

        std::string print_as_list();

        void canonicalize(){
            w.canonicalize();
            energy.canonicalize();
        }

        friend std::ostream& operator<<(std::ostream& os, const ScoreVector& score);
        friend bool operator==(const ScoreVector& a, const ScoreVector& b);
        friend bool operator!=(const ScoreVector& a, const ScoreVector& b);
    };

    mpq_class get_mpq_from_word(std::string word);

    class RNASequence {
        /**
           Representation of an RNA sequence
         **/
    public:
        RNASequence() {}; // Default constructor to satisfy compiler
        RNASequence(const std::string& seq); // Construct from a string encoding the sequence
        RNASequence(const fs::path& filename); // Construct from a FASTA file

        const int len() const; // Return the length of the sequence
        const int base(int i) const; // Return the base at position i of the sequence
        const std::string subsequence(int i, int j) const; // Return the subsequence starting at position i and ending at j
        const bool can_pair(int i, int j) const; // Return true if the bases at i and j are a valid pair

        const char& operator[](const int index) const; // Retrieve a single base using index notation
        friend std::ostream& operator<<(std::ostream& out, const RNASequence& sequence); // Output the sequence to an ostream

    protected:
        std::string seq_txt;
        void sanitize_string();
    };

    class RNASequenceWithTables: public RNASequence {
        /**
           Intermediate data store for NNTM dynamic programming results
        **/
    public:
        RNASequenceWithTables() {}; // Default constructor for compiler
        RNASequenceWithTables(const RNASequence& seq, mpq_class infinity_value = mpq_class(9999999999999)); // TODO: This is a really gross way to handle this! Find a better one.

        std::vector<mpq_class> W;
        std::vector< std::vector<mpq_class> > V;
        std::vector< std::vector<mpq_class> > VBI;
        std::vector< std::vector<mpq_class> > VM;
        std::vector< std::vector<mpq_class> > WM;
        std::vector< std::vector<mpq_class> > WMPrime;

        void print_debug();
    };

    class RNAStructure {
        /**
           Representation of an RNA secondary structure
         **/
    public:
        RNASequence seq;

        RNAStructure() {}; // Default constructor for compiler
        RNAStructure(const RNASequence& seq); // Construct a (blank) structure from a given sequence

        const int len() const; // Return the length of the sequence

        void mark_pair(int i, int j); // Record that (i, j) are paired
        void mark_d5(int i); // Record that i dangles from the 5' end of i+1
        void mark_d3(int i); // Record that i dangles from the 3' end of i-1

        const bool does_d5(int i) const; // Return true if i dangles from the 5' end of i+1
        const bool does_d3(int i) const; // Return true if i dangles from the 3' end of i-1

        const std::deque< std::pair<int, int> > pairs() const; // Return a deque of the pairs in the structure

        const std::string& string() const; // Return the structure as a string

        const char& operator[](const int index) const; // Retrieve a single base using index notation
        friend std::ostream& operator<<(std::ostream& out, const RNAStructure& structure); // Output this structure as an ostream

    protected:
        std::string structure_as_chars;
    };

    class RNAStructureWithScore: public RNAStructure {
        /**
Representation of an RNA secondary structure that has been assigned a score
        **/
    public:
        ScoreVector score;

        RNAStructureWithScore() {}; // Default constructor for compiler
        RNAStructureWithScore(const RNAStructure& structure, const ScoreVector& score);
        friend std::ostream& operator<<(std::ostream& out, const RNAStructureWithScore& structure); // Output this structure and its scores as an ostream
    };

    class RNAStructureTree: public RNAStructure {
    public:
        IntervalTreeNode root;

        RNAStructureTree() {}; // Default constructor for compiler
        RNAStructureTree(const RNAStructure& structure);
    };

    enum dangle_mode {
        NO_DANGLE = 0,
        CHOOSE_DANGLE = 1,
        BOTH_DANGLE = 2,
    };

    dangle_mode convert_to_dangle_mode(int n);

    enum RNA_base {
        BASE_A = 0,
        BASE_C = 1,
        BASE_G = 2,
        BASE_U = 3,
    };
}
#endif
