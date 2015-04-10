#ifndef _NNTM_CLASS_H_
#define _NNTM_CLASS_H_

#include "nndb_constants.h"
#include "pmfe_types.h"
#include <gmpxx.h>
#include <vector>

namespace pmfe{
    class NNTM {
    public:
        const NNDBConstants constants;
        const dangle_mode dangles;

        NNTM(NNDBConstants constants, dangle_mode dangles);

        RNASequenceWithTables energy_tables(const RNASequence& seq) const;
        mpq_class minimum_energy(const RNASequenceWithTables& seq) const;
        RNAStructureWithScore mfe_structure(const RNASequenceWithTables& seq) const;

    protected:
        mpq_class Ed3(int i, int j, int k, const RNASequenceWithTables& seq) const;
        mpq_class Ed5(int i, int j, int k, const RNASequenceWithTables& seq) const;
        mpq_class auPenalty(int i, int j, const RNASequenceWithTables& seq) const;
        mpq_class eL(int i, int j, int ip, int jp, const RNASequenceWithTables& seq) const;
        mpq_class eH(int i, int j, const RNASequenceWithTables& seq) const;
        mpq_class eS(int i, int j, const RNASequenceWithTables& seq) const;
        mpq_class calcVBI(int i, int j, const RNASequenceWithTables& seq) const;

        void traceW(int i, const RNASequenceWithTables& seq, RNAStructure& structure, ScoreVector& score) const;
        mpq_class traceV(int i, int j, const RNASequenceWithTables& seq, RNAStructure& structure, ScoreVector& score) const;
        mpq_class traceVM(int i, int j, const RNASequenceWithTables& seq, RNAStructure& structure, ScoreVector& score) const;
        mpq_class traceVBI(int i, int j, const RNASequenceWithTables& seq, RNAStructure& structure, ScoreVector& score) const;
        mpq_class traceWM(int i, int j, const RNASequenceWithTables& seq, RNAStructure& structure, ScoreVector& score) const;
        mpq_class traceWMPrime(int i, int j, const RNASequenceWithTables& seq, RNAStructure& structure, ScoreVector& score) const;

        const static int MAXLOOP = 30; /* The maximum loop size. */
        const static int TURN = 3; /* Minimum size of a hairpin loop. */
    };
};
#endif
