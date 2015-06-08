// Copyright (c) 2014 Andrew Gainer-Dewar.

#ifndef PARAMETRIZER_TYPES_H
#define PARAMETRIZER_TYPES_H

#include <string>
#include <utility>
#include <iostream>

#include <deque>
#include <stack>

#include <boost/filesystem.hpp>
#include "boost/multi_array.hpp"

#include "interval_tree.h"
#include "rational.h"

namespace pmfe {
    extern const Rational multiloop_default;
    extern const Rational unpaired_default;
    extern const Rational branch_default;
    extern const Rational dummy_default;

    namespace fs = boost::filesystem;

    enum RNA_base {
        BASE_A = 0,
        BASE_C = 1,
        BASE_G = 2,
        BASE_U = 3,
    };

    enum subopt_label {
        lW,
        lV,
        lVBI,
        lM,
        lM1,
    };

    enum dangle_mode {
        NO_DANGLE = 0,
        CHOOSE_DANGLE = 1,
        BOTH_DANGLE = 2,
    };

    class ParameterVector {
    public:
    ParameterVector(Rational multiloop_penalty = multiloop_default, Rational unpaired_penalty = unpaired_default, Rational branch_penalty = branch_default, Rational dummy_scaling = dummy_default) :
        multiloop_penalty(multiloop_penalty),
            unpaired_penalty(unpaired_penalty),
            branch_penalty(branch_penalty),
            dummy_scaling(dummy_scaling)
            {
            this->canonicalize();
        };
        Rational multiloop_penalty, unpaired_penalty, branch_penalty, dummy_scaling;


        void canonicalize() {
            multiloop_penalty.canonicalize();
            unpaired_penalty.canonicalize();
            branch_penalty.canonicalize();
            dummy_scaling.canonicalize();
        }

        std::string print_as_list();

        friend std::ostream& operator<<(std::ostream& os, const ParameterVector& params);
        friend bool operator==(const ParameterVector& a, const ParameterVector& b);
        friend bool operator!=(const ParameterVector& a, const ParameterVector& b);
    };

    class ScoreVector {
    public:
    ScoreVector(Integer multiloops = 0, Integer unpaired = 0, Integer branches = 0, Rational w = Rational(0), Rational energy = Rational(0)) :
        multiloops(multiloops),
            branches(branches),
            unpaired(unpaired),
            w(w),
            energy(energy)
            {
            this->canonicalize();
        };
        Integer multiloops, branches, unpaired;
        Rational w, energy;

        std::string print_as_list();

        void canonicalize(){
            w.canonicalize();
            energy.canonicalize();
        }

        friend std::ostream& operator<<(std::ostream& os, const ScoreVector& score);
        friend bool operator==(const ScoreVector& a, const ScoreVector& b);
        friend bool operator!=(const ScoreVector& a, const ScoreVector& b);
        friend bool operator<(const ScoreVector& a, const ScoreVector& b);

        ScoreVector& operator+=(const ScoreVector& rhs);
        friend ScoreVector operator+(const ScoreVector& lhs, const ScoreVector& rhs);
    };

    Rational get_rational_from_word(std::string word);

    class RNASequence {
        /**
           Representation of an RNA sequence
         **/
    public:
        RNASequence() {}; // Default constructor to satisfy compiler
        RNASequence(const std::string& seq); // Construct from a string encoding the sequence
        RNASequence(const fs::path& filename); // Construct from a FASTA file

        int len() const; // Return the length of the sequence
        int base(int i) const; // Return the base at position i of the sequence
        std::string subsequence(int i, int j) const; // Return the subsequence starting at position i and ending at j
        bool can_pair(int i, int j) const; // Return true if the bases at i and j are a valid pair

        char operator[](const int index) const; // Retrieve a single base using index notation
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
        RNASequenceWithTables(const RNASequence& seq);

        boost::multi_array<Rational, 1> W;
        boost::multi_array<Rational, 2> V;
        boost::multi_array<Rational, 2> VBI;
        boost::multi_array<Rational, 2> VM;
        boost::multi_array<Rational, 2> WM;
        boost::multi_array<Rational, 2> WMPrime;
        boost::multi_array<Rational, 2> FM;
        boost::multi_array<Rational, 2> FM1;

        bool energy_tables_populated = false;
        bool subopt_tables_populated = false;

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
        RNAStructure(const RNASequence& seq, const std::string& structure); // Construct a structure from dots-and-braces notation

        int len() const; // Return the length of the sequence

        void mark_pair(int i, int j); // Record that (i, j) are paired
        void mark_d5(int i); // Record that i dangles from the 5' end of i+1
        void mark_d3(int i); // Record that i dangles from the 3' end of i-1

        bool does_d5(int i) const; // Return true if i dangles from the 5' end of i+1
        bool does_d3(int i) const; // Return true if i dangles from the 3' end of i-1

        std::deque< std::pair<int, int> > pairs() const; // Return a deque of the pairs in the structure

        std::string string() const; // Return the structure as a new-style dots-and-braces string
        std::string old_string() const; // Return the structure as an old-style dots-and-braces string

        char operator[](const int index) const; // Retrieve a single base using index notation
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

        friend bool operator<(const RNAStructureWithScore& lhs, const RNAStructureWithScore& rhs);
    };

    class RNAStructureTree: public RNAStructure {
    public:
        IntervalTreeNode root;

        RNAStructureTree() {}; // Default constructor for compiler
        RNAStructureTree(const RNAStructure& structure);
    };

    class Segment {
        /**
           Representation of a segment in a suboptimal structure processing stack
         **/
    public:
        int i, j;
        subopt_label label;
        Rational minimum_energy;

    Segment(int i, int j, subopt_label label, Rational minimum_energy):
        i(i),
            j(j),
            label(label),
            minimum_energy(minimum_energy)
            {};

        friend std::ostream& operator<<(std::ostream& out, const Segment& seg); // Output this segment as an ostream
    };

    class RNAPartialStructure: public RNAStructure {
        /**
           Representation of a partial RNA secondary structure
        **/
    public:
        RNAPartialStructure(); // Default constructor for compiler
        RNAPartialStructure(const RNASequence& seq, Rational known_energy = 0); // Construct a (blank) structure from a given sequence with specified energy

        void accumulate(Rational energy); // Add to the known energy
        Rational total() const; // Return the known energy
        void push(const Segment& seg); // Push a segment onto the stack
        void pop(); // Remove a segment from the stack
        Segment top() const; // Retrieve the top segment of the stack
        bool empty() const; // True if the stack is empty

    protected:
        std::stack<Segment> seg_stack;
        Rational known_energy;
    };

    typedef std::stack<RNAPartialStructure> PartialStructureStack;

    dangle_mode convert_to_dangle_mode(int n);
}
#endif
