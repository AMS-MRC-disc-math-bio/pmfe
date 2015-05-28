#ifndef _NNTM_CLASS_H_
#define _NNTM_CLASS_H_

#include "nndb_constants.h"
#include "pmfe_types.h"
#include <gmpxx.h>
#include <vector>
#include "thread_pool.h"

namespace pmfe{
    class NNTM {
    public:
        const NNDBConstants& constants;
        const dangle_mode dangles;

        NNTM(const NNDBConstants& constants, dangle_mode dangles, SimpleThreadPool& thread_pool);

        RNASequenceWithTables energy_tables(const RNASequence& seq) const;
        mpq_class minimum_energy(RNASequenceWithTables& seq) const;
        RNAStructureWithScore mfe_structure(const RNASequenceWithTables& seq) const;

        ScoreVector score(const RNAStructure& structure, bool compute_w = true) const;

        std::vector<RNAStructureWithScore> suboptimal_structures(RNASequenceWithTables& seq, mpq_class delta, bool sorted = false) const;

    protected:
        // MFE helpers
        void populate_energy_tables(RNASequenceWithTables& seq) const;
        void populate_energy_tables(int i, int j, RNASequenceWithTables& seq) const;
        void populate_subopt_tables(RNASequenceWithTables& seq) const;
        void populate_subopt_tables(int i, int j, RNASequenceWithTables& seq) const;
        mpq_class Ed3(int i, int j, const RNASequence& seq, bool inside = false) const;
        mpq_class Ed5(int i, int j, const RNASequence& seq, bool inside = false) const;
        mpq_class auPenalty(int i, int j, const RNASequence& seq) const;
        mpq_class eLL(int size) const;
        mpq_class eL(int i, int j, int ip, int jp, const RNASequence& seq) const;
        mpq_class eH(int i, int j, const RNASequence& seq) const;
        mpq_class eS(int i, int j, const RNASequence& seq) const;
        mpq_class calcVBI(int i, int j, const RNASequenceWithTables& seq) const;

        // Traceback helpers
        void traceW(int i, const RNASequenceWithTables& seq, RNAStructure& structure, ScoreVector& score) const;
        mpq_class traceV(int i, int j, const RNASequenceWithTables& seq, RNAStructure& structure, ScoreVector& score) const;
        mpq_class traceVM(int i, int j, const RNASequenceWithTables& seq, RNAStructure& structure, ScoreVector& score) const;
        mpq_class traceVBI(int i, int j, const RNASequenceWithTables& seq, RNAStructure& structure, ScoreVector& score) const;
        mpq_class traceWM(int i, int j, const RNASequenceWithTables& seq, RNAStructure& structure, ScoreVector& score) const;
        mpq_class traceWMPrime(int i, int j, const RNASequenceWithTables& seq, RNAStructure& structure, ScoreVector& score) const;

        // Scoring helpers
        ScoreVector scoreTree(const RNAStructureTree& tree) const; // Score a whole structure tree
        ScoreVector scoreInternalNodeRecursively(const RNAStructureTree& tree, const IntervalTreeNode& node) const; // Score a given internal node of a structure tree, descending to its children
        ScoreVector scoreMUnpairedRegion(const RNAStructure& structure, int i1, int i2, int j1, int j2, bool is_external = false) const; // Compute the energy associated to an unpaired region in a multiloop or external loop
        ScoreVector scoreM(const RNAStructureTree& tree, const IntervalTreeNode& node) const; // Compute the energy associated to a multiloop node
        ScoreVector scoreE(const RNAStructureTree& tree) const; // Compute the energy associated to the external loop node

        // Suboptimal structure helpers
        bool subopt_process_top_structure(const RNASequenceWithTables& seq, RNAPartialStructure& ps, PartialStructureStack& pstack, mpq_class upper_bound) const;
        bool subopt_traceV(int i, int j, const RNASequenceWithTables& seq, RNAPartialStructure& ps, PartialStructureStack& pstack, mpq_class upper_bound) const;
        bool subopt_traceVBI(int i, int j, const RNASequenceWithTables& seq, RNAPartialStructure& ps, PartialStructureStack& pstack, mpq_class upper_bound) const;
        bool subopt_traceW(int i, int j, const RNASequenceWithTables& seq, RNAPartialStructure& ps, PartialStructureStack& pstack, mpq_class upper_bound) const;
        bool subopt_traceM1(int i, int j, const RNASequenceWithTables& seq, RNAPartialStructure& ps, PartialStructureStack& pstack, mpq_class upper_bound) const;
        bool subopt_traceM(int i, int j, const RNASequenceWithTables& seq, RNAPartialStructure& ps, PartialStructureStack& pstack, mpq_class upper_bound) const;

        // Threading
        SimpleThreadPool& thread_pool;

        // Configurable constants
        const static int MAXLOOP = 30; /* The maximum loop size. */
        const static int TURN = 3; /* Minimum size of a hairpin loop. */
    };
};
#endif
