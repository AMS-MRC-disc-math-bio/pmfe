#ifndef _NNTM_CLASS_H_
#define _NNTM_CLASS_H_

#include "nndb_constants.h"
#include "pmfe_types.h"
#include "rational.h"

#include <vector>


namespace pmfe{
    class NNTM {
    public:
        const NNDBConstants& constants;
        const dangle_mode dangles;

        NNTM(const NNDBConstants& constants, dangle_mode dangles);

        RNASequenceWithTables energy_tables(const RNASequence& seq) const;
        Rational minimum_energy(RNASequenceWithTables& seq) const;
        RNAStructureWithScore mfe_structure(const RNASequenceWithTables& seq) const;

        ScoreVector score(const RNAStructure& structure, bool compute_w = true) const;

        std::vector<RNAStructureWithScore> suboptimal_structures(RNASequenceWithTables& seq, Rational delta, bool sorted = false) const;

    protected:
        // MFE helpers
        void populate_energy_tables(RNASequenceWithTables& seq) const;
        void populate_energy_tables(int i, int j, RNASequenceWithTables& seq) const;
        void populate_subopt_tables(RNASequenceWithTables& seq) const;
        void populate_subopt_tables(int i, int j, RNASequenceWithTables& seq) const;
        Rational Ed3(int i, int j, const RNASequence& seq, bool inside = false) const;
        Rational Ed5(int i, int j, const RNASequence& seq, bool inside = false) const;
        Rational auPenalty(int i, int j, const RNASequence& seq) const;
        Rational eLL(int size) const;
        Rational eL(int i, int j, int ip, int jp, const RNASequence& seq) const;
        Rational eH(int i, int j, const RNASequence& seq) const;
        Rational eS(int i, int j, const RNASequence& seq) const;
        Rational calcVBI(int i, int j, const RNASequenceWithTables& seq) const;

        // Traceback helpers
        bool traceW(int i, const RNASequenceWithTables& seq, RNAStructure& structure, ScoreVector& score) const;
        Rational traceV(int i, int j, const RNASequenceWithTables& seq, RNAStructure& structure, ScoreVector& score) const;
        Rational traceVM(int i, int j, const RNASequenceWithTables& seq, RNAStructure& structure, ScoreVector& score) const;
        Rational traceVBI(int i, int j, const RNASequenceWithTables& seq, RNAStructure& structure, ScoreVector& score) const;
        Rational traceWM(int i, int j, const RNASequenceWithTables& seq, RNAStructure& structure, ScoreVector& score) const;
        Rational traceWMPrime(int i, int j, const RNASequenceWithTables& seq, RNAStructure& structure, ScoreVector& score) const;

        // Scoring helpers
        ScoreVector scoreTree(const RNAStructureTree& tree) const; // Score a whole structure tree
        ScoreVector scoreInternalNodeRecursively(const RNAStructureTree& tree, const IntervalTreeNode& node) const; // Score a given internal node of a structure tree, descending to its children
        ScoreVector scoreMUnpairedRegion(const RNAStructure& structure, int i1, int i2, int j1, int j2, bool is_external = false) const; // Compute the energy associated to an unpaired region in a multiloop or external loop
        ScoreVector scoreM(const RNAStructureTree& tree, const IntervalTreeNode& node) const; // Compute the energy associated to a multiloop node
        ScoreVector scoreE(const RNAStructureTree& tree) const; // Compute the energy associated to the external loop node

        // Suboptimal structure helpers
        bool subopt_process_top_structure(const RNASequenceWithTables& seq, RNAPartialStructure& ps, PartialStructureStack& pstack, Rational upper_bound) const;
        bool subopt_traceV(int i, int j, const RNASequenceWithTables& seq, RNAPartialStructure& ps, PartialStructureStack& pstack, Rational upper_bound) const;
        bool subopt_traceVBI(int i, int j, const RNASequenceWithTables& seq, RNAPartialStructure& ps, PartialStructureStack& pstack, Rational upper_bound) const;
        bool subopt_traceW(int i, int j, const RNASequenceWithTables& seq, RNAPartialStructure& ps, PartialStructureStack& pstack, Rational upper_bound) const;
        bool subopt_traceM1(int i, int j, const RNASequenceWithTables& seq, RNAPartialStructure& ps, PartialStructureStack& pstack, Rational upper_bound) const;
        bool subopt_traceM(int i, int j, const RNASequenceWithTables& seq, RNAPartialStructure& ps, PartialStructureStack& pstack, Rational upper_bound) const;

        // Configurable constants
        const static int MAXLOOP = 30; /* The maximum loop size. */
        const static int TURN = 3; /* Minimum size of a hairpin loop. */
    };
};
#endif
