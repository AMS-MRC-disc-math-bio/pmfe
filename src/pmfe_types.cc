// Copyright (c) 2014 Andrew Gainer-Dewar.

#include <utility>
#include <vector>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <algorithm>

#include <set>
#include <deque>

#include <assert.h>

#include "pmfe_types.h"
#include "rational.h"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>

#include <boost/assign/list_of.hpp>

namespace pmfe {
    namespace fs = boost::filesystem;

    const Rational multiloop_default = Rational(17, 5);
    const Rational unpaired_default = Rational(0);
    const Rational branch_default = Rational(2, 5);
    const Rational dummy_default = Rational(1);

    const char d3symb = '>';
    const char d5symb = '<';
    const char p3symb = ')';
    const char p5symb = '(';
    const char blanksymb = '.';

    std::ostream& operator<<(std::ostream& os, const ParameterVector& params) {
        os << "Multiloop penalty: " << params.multiloop_penalty.get_str(10) << " ≈ " << params.multiloop_penalty.get_d() << std::endl
           << "Unpaired penalty: " << params.unpaired_penalty.get_str(10) << " ≈ " << params.unpaired_penalty.get_d() << std::endl
           << "Branch penalty: " << params.branch_penalty.get_str(10) << " ≈ " << params.branch_penalty.get_d() << std::endl
           << "Dummy scaling: " << params.dummy_scaling.get_str(10) << " ≈ " << params.dummy_scaling.get_d() << std::endl;
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

    std::string ParameterVector::print_as_list() {
        std::ostringstream res;
        res << "["
            << multiloop_penalty.get_str(10)
            << ", "
            << unpaired_penalty.get_str(10)
            << ", "
            << branch_penalty.get_str(10)
            << ", "
            << dummy_scaling.get_str(10)
            << "]";
        return res.str();
    };

    std::ostream& operator<<(std::ostream& os, const ScoreVector& score) {
        os << "Multiloops: " << score.multiloops.get_str(10) << std::endl
           << "Unpaired bases: " << score.unpaired.get_str(10) << std::endl
           << "Branches: " << score.branches.get_str(10) << std::endl
           << "w: " << score.w.get_str(10) << " ≈ " << score.w.get_d() << std::endl
           << "Parametrized energy: " << score.energy.get_str(10) << " ≈ " << score.energy.get_d() << std::endl;
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

    bool operator<(const ScoreVector& a, const ScoreVector& b) {
        // The most useful ordering on score vectors is lexicographic with energy first
        std::vector<Rational> a_vect = {a.energy, a.multiloops, a.unpaired, a.branches, a.w};
        std::vector<Rational> b_vect = {b.energy, b.multiloops, b.unpaired, b.branches, b.w};

        return std::lexicographical_compare(a_vect.begin(), a_vect.end(), b_vect.begin(), b_vect.end());
    }

    ScoreVector& ScoreVector::operator+=(const ScoreVector& rhs) {
        this->multiloops += rhs.multiloops;
        this->unpaired += rhs.unpaired;
        this->branches += rhs.branches;
        this->w += rhs.branches;
        this->energy += rhs.energy;

        return *this;
    }

    ScoreVector operator+(const ScoreVector& lhs, const ScoreVector& rhs) {
        ScoreVector result = lhs;
        result += rhs;
        return result;
    }

    std::string ScoreVector::print_as_list() {
        std::ostringstream res;
        res << "["
            << multiloops.get_str(10)
            << ", "
            << unpaired.get_str(10)
            << ", "
            << branches.get_str(10)
            << ", "
            << w.get_str(10)
            << "]";
        return res.str();
    };

    Rational get_rational_from_word(std::string word) {
        Rational result;
        std::size_t decimalpoint = word.find('.');
        bool negative = (word.find('-') != std::string::npos);
        if (decimalpoint != std::string::npos) {
            std::string intpart = word.substr(0, decimalpoint);
            std::string fracpart = word.substr(decimalpoint+1);

            if (intpart == "") intpart = "0";
            Integer theint (intpart, 10);
            theint = abs(theint);

            // Carve out the fractional part. Surprisingly fiddly!
            int fracdenom = pow(10, fracpart.length());
            Integer fracval (fracpart, 10);
            Rational thefrac (fracval, fracdenom);
            thefrac.canonicalize();

            result = theint + thefrac;
            if (negative) result *= -1;
        } else {
            Rational thevalue (word);
            result = thevalue;
        }
        result.canonicalize();
        return result;
    };

    RNASequence::RNASequence(const std::string& seq):
        seq_txt(seq)
    {
        valid_pairs.resize(boost::extents[len()][len()]);
        preprocess();
    };

    RNASequence::RNASequence(const fs::path& filename) {
        if (not fs::is_regular_file(filename)) {
            std::stringstream error_message;
            error_message << "Path " << filename << " does not point to a valid file." ;
            throw std::invalid_argument(error_message.str());
        }

        std::string tempseq = "";
        fs::ifstream filestream (filename);

        std::string line;
        while (std::getline(filestream, line)) {
            // Remove any DOS-style line terminators
            line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

            // As well as any spaces or tabs
            line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end()); // We use the C implementation of isspace to avoid goofy template problems

            // exclude lines starting with FASTA comment characters
            if(line[0] != ';' && line[0] != '>' && line.length() > 0)
                tempseq += line;
        }

        seq_txt = tempseq;

        valid_pairs.resize(boost::extents[len()][len()]);
        preprocess();
    }

    void RNASequence::preprocess() {
        // Apply some processing to each character of the seq_txt string
        for (int i = 0; i < len(); ++i) {
            // Capitalize if the base is valid, throw an exception if not
            switch(base(i)) {
            case BASE_A:
                seq_txt[i] = 'A';
                break;

            case BASE_C:
                seq_txt[i] = 'C';
                break;

            case BASE_G:
                seq_txt[i] = 'G';
                break;

            case BASE_U:
                seq_txt[i] = 'U';
                break;

            default:
                std::stringstream error_message;
                error_message << "At position " << i << ", found " << seq_txt[i] << ", which is an invalid RNA base.";
                throw std::invalid_argument(error_message.str());
                break;
            }
        }

        // Populate valid_pairs
        std::set<std::string> RNAPairs = {"AU", "UA", "CG", "GC", "GU", "UG"};

        for (int i = 0; i < len(); ++i) {
            for (int j = 0; j < i; ++j) {
                std::string thepair = seq_txt.substr(i, 1) + seq_txt.substr(j, 1);
                if (RNAPairs.count(thepair) == 0) {
                    valid_pairs[i][j] = valid_pairs[j][i] = false;
                } else {
                    valid_pairs[i][j] = valid_pairs[j][i] = true;
                }
            }
        }
    }

    int RNASequence::len() const {
        return seq_txt.length();
    }

    int RNASequence::base(int i) const {
        char base = seq_txt[i];
        switch(base) {
        case 'A':
        case 'a':
            return BASE_A;
        case 'C':
        case 'c':
            return BASE_C;
        case 'G':
        case 'g':
            return BASE_G;
        case 'U':
        case 'u':
        case 'T':
        case 't':
            return BASE_U;
        default:
            return -1;
        }
    }

    std::string RNASequence::subsequence(int i, int j) const {
        // Return the subsequence starting at i and ending at j (inclusive)
        assert(i <= j);
        return seq_txt.substr(i, j-i+1);
    }

    bool RNASequence::can_pair(int i, int j) const {
        return valid_pairs[i][j];
    }

    char RNASequence::operator[](const int index) const {
        return seq_txt[index];
    }

    std::ostream& operator<<(std::ostream& os, const RNASequence& sequence) {
        os << sequence.seq_txt;
        return os;
    }

    RNAStructure::RNAStructure(const RNASequence& seq):
        seq(seq) {
        structure_as_chars = std::string(len(), blanksymb);
    };

    RNAStructure::RNAStructure(const RNASequence& seq, const std::string& structure):
        // TODO: Sanity-check structure!
        seq(seq),
        structure_as_chars(structure)
    {};

    char RNAStructure::operator[](const int index) const {
        return seq[index];
    }

    std::string RNAStructure::string() const {
        return structure_as_chars;
    }

    std::string RNAStructure::old_string() const {
        std::string result = structure_as_chars;
        std::replace(result.begin(), result.end(), d3symb, blanksymb);
        std::replace(result.begin(), result.end(), d5symb, blanksymb);
        return result;
    }

    RNASequenceWithTables::RNASequenceWithTables(const RNASequence& seq):
        RNASequence(seq),
        W(boost::extents[len()]),
        V(boost::extents[len()][len()]),
        VBI(boost::extents[len()][len()]),
        VM(boost::extents[len()][len()]),
        WM(boost::extents[len()][len()]),
        WMPrime(boost::extents[len()][len()]),
        FM(boost::extents[len()][len()]),
        FM1(boost::extents[len()][len()])
    {
        // All the arrays should be filled with infinity
        std::fill(W.data(), W.data() + W.num_elements(), Rational::infinity());
        std::fill(V.data(), V.data() + V.num_elements(), Rational::infinity());
        std::fill(VBI.data(), VBI.data() + VBI.num_elements(), Rational::infinity());
        std::fill(VM.data(), VM.data() + VM.num_elements(), Rational::infinity());
        std::fill(WM.data(), WM.data() + WM.num_elements(), Rational::infinity());
        std::fill(WMPrime.data(), WMPrime.data() + WMPrime.num_elements(), Rational::infinity());
        std::fill(FM.data(), FM.data() + FM.num_elements(), Rational::infinity());
        std::fill(FM1.data(), FM1.data() + FM1.num_elements(), Rational::infinity());
    }

    void RNASequenceWithTables::print_debug(){
        printf("Intermediate tables:\n");
        for (int b = 4; b <= len(); ++b) {
            for (int i = 1; i <= len() - b; ++i) {
                int j = i + b;
                printf("VBI[%d][%d] = %f\n", i, j, VBI[i][j].get_d());
                printf("VM[%d][%d] = %f\n", i, j, VM[i][j].get_d());
                printf("V[%d][%d] = %f\n", i, j, V[i][j].get_d());
                printf("WMPrime[%d][%d] = %f\n", i, j, WMPrime[i][j].get_d());
                printf("WM[%d][%d] = %f\n", i, j, WM[i][j].get_d());
            }
        }

        printf("\nW final W table:\n");
        for (int i = 1; i <= len(); ++i) {
            printf("W[%d] = %f\n", i, W[i].get_d());
        }
    }

    int RNAStructure::len() const {
        return seq.len();
    }

    void RNAStructure::mark_pair(int i, int j) {
        assert (structure_as_chars[i] == blanksymb);
        assert (structure_as_chars[j] == blanksymb);

        structure_as_chars[i] = p5symb;
        structure_as_chars[j] = p3symb;
    }

    void RNAStructure::mark_d5(int i) {
        assert (structure_as_chars[i] == blanksymb);
        structure_as_chars[i] = d5symb;
    }

    void RNAStructure::mark_d3(int i) {
        assert (structure_as_chars[i] == blanksymb);
        structure_as_chars[i] = d3symb;
    }

    bool RNAStructure::does_d5(int i) const {
        return (structure_as_chars[i] == d5symb);
    }

    bool RNAStructure::does_d3(int i) const {
        return (structure_as_chars[i] == d3symb);
    }

    std::deque< std::pair<int, int> > RNAStructure::pairs() const {
        /*
          Return the pairs in this structure, sorted in increasing order of first base
        */

        // We find the pairs using a simple search algorithm.
        // We read through the structure string once.
        // Each time we find an open parenthesis, we add its index to the deque 'starts'.
        // Each time we find a close parenthesis, we pop the last start off the deque,
        // form a pair, and add it to the results.

        std::deque<int> starts;
        std::deque< std::pair<int, int> > results;
        for (unsigned int i = 0; i < structure_as_chars.length(); ++i) {
            char c = structure_as_chars[i];
            switch (c) {
            case p5symb:
            {
                starts.push_back(i);
                break;
            }

            case p3symb:
            {
                int start = starts.back();
                starts.pop_back();
                results.push_back(std::make_pair(start, i));
                break;
            }

            default:
                break;
            }
        }

        std::sort(results.begin(), results.end());
        return results;
    }

    std::ostream& operator<<(std::ostream& os, const RNAStructure& structure) {
        os << structure.structure_as_chars;
        return os;
    }

    RNAStructureWithScore::RNAStructureWithScore(const RNAStructure& structure, const ScoreVector& score):
        RNAStructure(structure), score(score) {};

    std::ostream& operator<<(std::ostream& os, const RNAStructureWithScore& structure) {
        os << structure.structure_as_chars << "\t"
           << structure.score.multiloops << "\t"
           << structure.score.unpaired << "\t"
           << structure.score.branches << "\t"
           << structure.score.w << "\t"
           << structure.score.energy;
        return os;
    }

    bool operator<(const RNAStructureWithScore& lhs, const RNAStructureWithScore& rhs) {
        return lhs.score < rhs.score;
    }

    RNAStructureTree::RNAStructureTree(const RNAStructure& structure):
        RNAStructure(structure)
    {
        root = IntervalTreeNode(-1, structure.len());
        std::deque< std::pair<int, int> > pairs = structure.pairs();
        for (std::deque< std::pair<int, int> >::const_iterator pair = pairs.begin(); pair != pairs.end(); ++pair) {
            root.insert(pair->first, pair->second);
        }
    };

    std::ostream& operator<<(std::ostream& os, const Segment& seg) {
        os << "(" << seg.i << ", " << seg.j << ")_" << seg.label;
        return os;
    }

    RNAPartialStructure::RNAPartialStructure():
        RNAStructure(),
        known_energy(0)
    {};

    RNAPartialStructure::RNAPartialStructure(const RNASequence& seq, Rational known_energy):
        RNAStructure(seq),
        known_energy(known_energy)
    {};

    void RNAPartialStructure::accumulate(Rational energy) {
        known_energy += energy;
    };

    Rational RNAPartialStructure::total() const {
        return known_energy;
    };

    void RNAPartialStructure::push(const Segment& seg) {
        known_energy += seg.minimum_energy;
        seg_stack.push(seg);
    };

    void RNAPartialStructure::pop() {
        known_energy -= top().minimum_energy;
        seg_stack.pop();
    };

    Segment RNAPartialStructure::top() const {
        return seg_stack.top();
    };

    bool RNAPartialStructure::empty() const {
        return seg_stack.empty();
    };

    dangle_mode convert_to_dangle_mode(int n) {
        switch (n) {
        case 0:
            return NO_DANGLE;
            break;
        case 1:
            return CHOOSE_DANGLE;
            break;
        case 2:
            return BOTH_DANGLE;
            break;
        default:
            exit(EXIT_FAILURE);
            break;
        }
    }
}
