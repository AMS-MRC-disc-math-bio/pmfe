// Copyright (c) 2015 Andrew Gainer-Dewar.

#include <stdlib.h>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <iostream>

#include "nntm.h"
#include "nndb_constants.h"
#include "pmfe_types.h"

#include <gmpxx.h>
#include "boost/multi_array.hpp"

namespace pmfe {
    void NNTM::populate_subopt_tables(RNASequenceWithTables& seq) const {
        /*
          Construct the helper tables for the subopt traceback algorithm
        */
        // Input specification
        assert(not seq.subopt_tables_populated);

        // Populate FM1 and FM
        for (int b = TURN+1; b <= seq.len() - 1; ++b) {
            SimpleJobGroup job_group(thread_pool);
            for (int i = 0; i <= seq.len() - 1 - b; ++i) {
                job_group.post(boost::bind(&NNTM::populate_subopt_tables, this, i, i+b, std::ref(seq)));
            }

            job_group.wait_for_all_jobs();
        }
        seq.subopt_tables_populated = true;
    }

    void NNTM::populate_subopt_tables(int i, int j, RNASequenceWithTables& seq) const {
        // FM begin
        std::deque<mpq_class> fm1_vals;
        fm1_vals.push_back(constants.INFINITY_);

        int minl = i+TURN+1;
        for (int l = minl; l <= j; ++l) {
            switch (dangles) {

            case NO_DANGLE:
            {
                fm1_vals.push_back(seq.V[i][i] + auPenalty(i, l, seq) + constants.multConst[1]*(j-l) + constants.multConst[2]);
                break;
            }

            case CHOOSE_DANGLE:
            {
                mpq_class d5 = Ed5(i+1, l, seq);
                mpq_class d3 = Ed3(i, l-1, seq);
                mpq_class d53 = Ed5(i+1, l-1, seq) + Ed3(i+1, l-1, seq);

                fm1_vals.push_back(seq.V[i][l] + auPenalty(i, l, seq) + constants.multConst[1] * (j-l) + constants.multConst[2]);

                if (l >= minl + 1)
                    fm1_vals.push_back(seq.V[i+1][l] + auPenalty(i+1, l, seq) + d5 + constants.multConst[1] * (j-l+1) + constants.multConst[2]);

                if (l >= minl + 1)
                    fm1_vals.push_back(seq.V[i][l-1] + auPenalty(i, l-1, seq) + d3 + constants.multConst[1] * (j-(l-1)) + constants.multConst[2]);

                if (l >= minl + 2)
                    fm1_vals.push_back(seq.V[i+1][l-1] + auPenalty(i+1, l-1, seq) + d53 + constants.multConst[1]* (j-(l-1)+1) + constants.multConst[2]);

                break;
            }

            case BOTH_DANGLE:
            {
                mpq_class d5 = Ed5(i, l, seq);
                mpq_class d3 = Ed3(i, l, seq);
                fm1_vals.push_back(seq.V[i][l] + auPenalty(i, l, seq) + d5 + d3 + constants.multConst[1] * (j-l) + constants.multConst[2]);
                break;
            }

            default:
            {
                exit(EXIT_FAILURE);
                break;
            }
            }
        }
        seq.FM1[i][j] = *std::min_element(fm1_vals.begin(), fm1_vals.end());

        std::deque<mpq_class> fm_vals;
        fm_vals.push_back(constants.INFINITY_);

        for (int k = i+TURN+1; k <= j-TURN-1; ++k) {
            fm_vals.push_back(seq.FM[i][k-1] + seq.FM1[k][j]);
        }

        for (int k = i; k <= j-TURN-1; ++k) {
            fm_vals.push_back(seq.FM1[k][j] + constants.multConst[1]*(k-i));
        }
        seq.FM[i][j] = *std::min_element(fm_vals.begin(), fm_vals.end());
    }

    std::vector<RNAStructureWithScore> NNTM::suboptimal_structures(RNASequenceWithTables& seq, mpq_class delta, bool sorted) const {
        // Ensure tables are available
        if (not seq.subopt_tables_populated) {
            populate_subopt_tables(seq);
        }

        // Set up variables
        mpq_class mfe = minimum_energy(seq);
        mpq_class upper_bound = mfe + delta;

        PartialStructureStack pstack;
        std::vector<RNAStructureWithScore> possible_structures;

        // Construct the initial partial sequence and add it to the stack
        RNAPartialStructure first(seq);
        first.push(Segment(0, seq.len()-1, lW, mfe));
        pstack.push(first);

        // Main processing loop
        while (!pstack.empty()) {
            RNAPartialStructure ps = pstack.top();
            pstack.pop();

            if (ps.empty() ) {
                // In this case, this structure is fully evaluated
                // Score the structure
                RNAStructure structure = ps;
                ScoreVector score = this->score(structure);
                RNAStructureWithScore result(structure, score);

                if (ps.total() != score.energy) {
                    throw std::logic_error("Inconsistent energy in suboptimal structure calculation.");
                }

                if (ps.total() > upper_bound) {
                    throw std::logic_error("Invalid energy in suboptimal structure calculation.");
                }

                possible_structures.push_back(result);
            } else {
                // Otherwise, we need to process the structure
                bool pushed_something = subopt_process_top_structure(seq, ps, pstack, upper_bound);

                // If nothing was pushed to the stack, we still need to consider the rest of the partial-structure stack
                if (!pushed_something) {
                    pstack.push(ps);
                }
            }
        }

        if (sorted) {
            std::sort(possible_structures.begin(), possible_structures.end());
        }

        return possible_structures;
    }

    bool NNTM::subopt_process_top_structure(const RNASequenceWithTables& seq, RNAPartialStructure& ps, PartialStructureStack& pstack, mpq_class upper_bound) const {
        // Take the top structure from the stack
        Segment seg = ps.top();
        ps.pop();

        // If this segment is too short to contain a substructure, exit immediately
        if (seg.j - seg.i <= TURN) {
            return false;
        }

        // Otherwise, apply the appropriate traceback function
        bool pushed_something;
        switch (seg.label){
        case lW:
            pushed_something = subopt_traceW(seg.i, seg.j, seq, ps, pstack, upper_bound);
            break;

        case lV:
            pushed_something = subopt_traceV(seg.i, seg.j, seq, ps, pstack, upper_bound);
            break;

        case lVBI:
            pushed_something = subopt_traceVBI(seg.i, seg.j, seq, ps, pstack, upper_bound);
            break;

        case lM:
            pushed_something = subopt_traceM(seg.i, seg.j, seq, ps, pstack, upper_bound);
            break;

        case lM1:
            pushed_something = subopt_traceM1(seg.i, seg.j, seq, ps, pstack, upper_bound);
            break;

        default:
            throw std::logic_error("Invalid label on suboptimal segment.");
            break;
        }
        return pushed_something;
    };

    bool NNTM::subopt_traceV(int i, int j, const RNASequenceWithTables& seq, RNAPartialStructure& ps, PartialStructureStack& pstack, mpq_class upper_bound) const {
        // Input specification
        assert (0 <= i);
        assert (i <= j);
        assert (j < ps.len());

        bool pushed_something = false;

        // Hairpin Loop
        if (eH(i, j, seq) + ps.total() <= upper_bound) {
            RNAPartialStructure new_ps(ps);
            new_ps.accumulate(eH(i, j, seq));
            new_ps.mark_pair(i, j);
            pstack.push(new_ps);
            mpq_class new_t = new_ps.total();
            pushed_something = true;
        }

        // Stack
        if (eS(i, j, seq) + seq.V[i+1][j-1] + ps.total() <= upper_bound) {
            RNAPartialStructure new_ps(ps);
            new_ps.push(Segment(i+1, j-1, lV, seq.V[i+1][j-1]));
            new_ps.accumulate(eS(i, j, seq));
            new_ps.mark_pair(i, j);
            pstack.push(new_ps);
            pushed_something = true;
        }

        // Internal Loop
        if (seq.VBI[i][j] + ps.total() <= upper_bound) {
            RNAPartialStructure new_ps(ps);
            new_ps.push(Segment(i, j, lVBI, seq.VBI[i][j]));
            pstack.push(new_ps);
            pushed_something = true;
        }

        // Multiloop
        for (int k = i + 2; k <= j-TURN-1; ++k) {
            switch (dangles) {
            case NO_DANGLE:
            {
                mpq_class kenergy1 = seq.FM[i+1][k] + seq.FM1[k+1][j-1];
                mpq_class kenergy2 = auPenalty(i, j, seq) + constants.multConst[0] + constants.multConst[2];
                mpq_class kenergy_total = kenergy1 + kenergy2;
                if (kenergy_total + ps.total() <= upper_bound) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(i+1, k, lM, seq.FM[i+1][k]));
                    new_ps.push(Segment(k+1, j-1, lM1, seq.FM1[k+1][j-1]));
                    new_ps.accumulate(kenergy2);
                    new_ps.mark_pair(i, j);
                    pstack.push(new_ps);
                    pushed_something = true;
                }
                break;
            }

            case CHOOSE_DANGLE:
                // In CHOOSE_DANGLE mode, we need to consider dangles on the initiating pair of a multiloop
            {
                mpq_class d5 = Ed5(i, j, seq, true);
                mpq_class d3 = Ed3(i, j, seq, true);
                mpq_class d53 = d5 + d3;
                if (seq.FM[i+1][k] + seq.FM1[k+1][j-1] + auPenalty(i, j, seq) + constants.multConst[0] + constants.multConst[2] + ps.total() <= upper_bound) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(i+1, k, lM, seq.FM[i+1][k]));
                    new_ps.push(Segment(k+1, j-1, lM1, seq.FM1[k+1][j-1]));
                    new_ps.accumulate(auPenalty(i, j, seq) + constants.multConst[0] + constants.multConst[2]);
                    new_ps.mark_pair(i, j);
                    pstack.push(new_ps);
                    pushed_something = true;
                }
                if (k > i+2 && seq.FM[i+2][k] + seq.FM1[k+1][j-1] + auPenalty(i, j, seq) + d5 + constants.multConst[0] + constants.multConst[1] + constants.multConst[2] + ps.total() <= upper_bound) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(i+2, k, lM, seq.FM[i+2][k]));
                    new_ps.push(Segment(k+1, j-1, lM1, seq.FM1[k+1][j-1]));
                    new_ps.accumulate(auPenalty(i, j, seq) + d5 + constants.multConst[0] + constants.multConst[1] + constants.multConst[2]);
                    new_ps.mark_pair(i, j);
                    new_ps.mark_d3(i+1);
                    pstack.push(new_ps);
                    pushed_something = true;
                }
                if (k <= j-TURN-2 && seq.FM[i+1][k] + seq.FM1[k+1][j-2] + auPenalty(i, j, seq) + d3 + constants.multConst[0] + constants.multConst[1] + constants.multConst[2] + ps.total() <= upper_bound) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(i+1, k, lM, seq.FM[i+1][k]));
                    new_ps.push(Segment(k+1, j-2, lM1, seq.FM1[k+1][j-2]));
                    new_ps.accumulate(auPenalty(i, j, seq) + d3 + constants.multConst[0] + constants.multConst[1] + constants.multConst[2]);
                    new_ps.mark_pair(i, j);
                    new_ps.mark_d5(j-1);
                    pstack.push(new_ps);
                    pushed_something = true;
                }
                if (k > i+2 && k <= j-TURN-2 && seq.FM[i+2][k] + seq.FM1[k+1][j-2] + auPenalty(i, j, seq) + d53 + constants.multConst[0] + 2*constants.multConst[1] + constants.multConst[2] + ps.total() <= upper_bound) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(i+2, k, lM, seq.FM[i+2][k]));
                    new_ps.push(Segment(k+1, j-2, lM1, seq.FM1[k+1][j-2]));
                    new_ps.accumulate(auPenalty(i, j, seq) + d53 + constants.multConst[0] + 2*constants.multConst[1] + constants.multConst[2]);
                    new_ps.mark_pair(i, j);
                    new_ps.mark_d3(i+1);
                    new_ps.mark_d5(j-1);
                    pstack.push(new_ps);
                    pushed_something = true;
                }
                break;
            }

            case BOTH_DANGLE:
            {
                mpq_class d5 = Ed5(i, j, seq, true);
                mpq_class d3 = Ed3(i, j, seq, true);
                mpq_class kenergy1 = seq.FM[i+1][k] + seq.FM1[k+1][j-1];
                mpq_class kenergy2 = d5 + d3 + auPenalty(i, j, seq) + constants.multConst[0] + constants.multConst[2];
                mpq_class kenergy_total = kenergy1 + kenergy2;
                if (kenergy_total + ps.total() <= upper_bound) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(i+1, k, lM, seq.FM[i+1][k]));
                    new_ps.push(Segment(k+1, j-1, lM1, seq.FM1[k+1][j-1]));
                    new_ps.accumulate(kenergy2);
                    new_ps.mark_pair(i, j);
                    pstack.push(new_ps);
                    pushed_something = true;
                }
                break;
            }

            default:
            {
                exit(EXIT_FAILURE);
                break;
            }
            };
        }

        return pushed_something;
    }

    bool NNTM::subopt_traceVBI(int i, int j, const RNASequenceWithTables& seq, RNAPartialStructure& ps, PartialStructureStack& pstack, mpq_class upper_bound) const {
        // Input specification
        assert (0 <= i);
        assert (i < j);
        assert (j < ps.len());

        bool pushed_something = false;
        for (int p = i+1; p <= std::min(j-2-TURN, i+MAXLOOP+1) ; ++p) {
            int minq = j-i+p-MAXLOOP-2;
            if (minq < p+1+TURN)
                minq = p+1+TURN;

            int maxq;
            if (p == i+1) {
                maxq = j-2;
            } else {
                maxq = j-1;
            }

            for (int q = minq; q <= maxq; ++q) {
                if (seq.V[p][q] + eL(i, j, p, q, seq) + ps.total() <= upper_bound) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(p, q, lV, seq.V[p][q]));
                    new_ps.mark_pair(i, j);
                    new_ps.accumulate(eL(i, j, p, q, seq));
                    pstack.push(new_ps);
                    pushed_something = true;
                }
            }
        }

        return pushed_something;
    }

    // Wuchty case E = F
    bool NNTM::subopt_traceW(int i, int j, const RNASequenceWithTables& seq, RNAPartialStructure& ps, PartialStructureStack& pstack, mpq_class upper_bound) const {
        // Input specification
        assert (i == 0);
        assert (i < j);
        assert (j < ps.len());

        bool pushed_something = false;
        for (int l = i; l < j-TURN; ++l) {
            mpq_class wim1;
            if (l > 0) {
                wim1 = seq.W[l-1];
            } else {
                wim1 = 0;
            }

            switch (dangles){
            case NO_DANGLE:
            {
                mpq_class bonus = auPenalty(l, j, seq);
                if (seq.V[l][j] + wim1 + bonus + ps.total() <= upper_bound ) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(l, j, lV, seq.V[l][j]));
                    if (l > i) new_ps.push(Segment(i, l-1, lW, wim1));
                    new_ps.accumulate(bonus);
                    pstack.push(new_ps);
                    pushed_something = true;
                }
                break;
            }

            case CHOOSE_DANGLE:
            {
                mpq_class d5 = Ed5(l+1, j, seq);
                mpq_class d3 = Ed3(l, j-1, seq);
                mpq_class d53 = Ed5(l+1, j-1, seq) + Ed3(l+1, j-1, seq);

                if (seq.V[l][j] + auPenalty(l, j, seq) + wim1 + ps.total() <= upper_bound) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(l, j, lV, seq.V[l][j]));
                    if (l > i) new_ps.push(Segment(i, l-1, lW, wim1));
                    new_ps.accumulate(auPenalty(l, j, seq));
                    pstack.push(new_ps);
                    pushed_something = true;
                }
                if (l+1 < j-TURN && seq.V[l+1][j] + auPenalty(l+1, j, seq) + d5 + wim1 + ps.total() <= upper_bound) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(l+1, j, lV, seq.V[l+1][j]));
                    new_ps.mark_d5(l);
                    if (l > i) new_ps.push(Segment(i, l-1, lW, wim1));
                    new_ps.accumulate(auPenalty(l+1, j, seq) + d5);
                    pstack.push(new_ps);
                    pushed_something = true;
                }
                if (l < j-TURN-1 && seq.V[l][j-1] + auPenalty(l, j-1, seq) +  d3 + wim1 + ps.total() <= upper_bound) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(l, j-1, lV, seq.V[l][j-1]));
                    new_ps.mark_d3(j);
                    if (l > i) new_ps.push(Segment(i, l-1, lW, wim1));
                    new_ps.accumulate(auPenalty(l, j-1, seq) + d3);
                    pstack.push(new_ps);
                    pushed_something = true;
                }
                if (l+1 < j-TURN-1 && seq.V[l+1][j-1] + auPenalty(l+1, j-1, seq) + d53 + wim1 + ps.total() <= upper_bound) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(l+1, j-1, lV, seq.V[l+1][j-1]));
                    new_ps.mark_d5(l);
                    new_ps.mark_d3(j);
                    if (l > i) new_ps.push(Segment(i, l-1, lW, wim1));
                    new_ps.accumulate(auPenalty(l+1, j-1, seq) + d53);
                    pstack.push(new_ps);
                    pushed_something = true;
                }
                break;
            }

            case BOTH_DANGLE:
            {
                mpq_class bonus = auPenalty(l, j, seq);

                if (l > i) {
                    bonus += Ed5(l, j, seq);
                }

                if (j < seq.len() - 1) {
                    bonus += Ed3(l, j, seq);
                }

                if (seq.V[l][j] + wim1 + bonus + ps.total() <= upper_bound ) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(l, j, lV, seq.V[l][j]));
                    if (l > i) new_ps.push(Segment(i, l-1, lW, wim1));
                    new_ps.accumulate(bonus);
                    pstack.push(new_ps);
                    pushed_something = true;
                }
                break;
            }

            default:
                exit(EXIT_FAILURE);
                break;
            };
        };

        if (seq.W[j-1] + ps.total() <= upper_bound) {
            RNAPartialStructure new_ps(ps);
            new_ps.push(Segment(i, j-1, lW, seq.W[j-1]));
            pstack.push(new_ps);
            pushed_something = true;
        }

        return pushed_something;
    }

    bool NNTM::subopt_traceM1(int i, int j, const RNASequenceWithTables& seq, RNAPartialStructure& ps, PartialStructureStack& pstack, mpq_class upper_bound) const {
        // Input specification
        assert (0 <= i);
        assert (i < j);
        assert (j < ps.len());

        bool pushed_something = false;

        if (seq.FM1[i][j-1] + constants.multConst[1] + ps.total() <= upper_bound) {
            RNAPartialStructure new_ps(ps);
            new_ps.push(Segment(i, j-1, lM1, seq.FM1[i][j-1]));
            new_ps.accumulate(constants.multConst[1]);
            pstack.push(new_ps);
            pushed_something = true;
        }

        switch (dangles) {
        case NO_DANGLE:
        {
            mpq_class bonus = auPenalty(i, j, seq) + constants.multConst[2];
            if (seq.V[i][j] + bonus + ps.total() <= upper_bound) {
                RNAPartialStructure new_ps(ps);
                new_ps.push(Segment(i, j, lV, seq.V[i][j]));
                new_ps.accumulate(bonus);
                pstack.push(new_ps);
                pushed_something = true;
            }
            break;
        }

        case CHOOSE_DANGLE:
        {
            mpq_class d5 = Ed5(i+1, j, seq);
            mpq_class d3 = Ed3(i, j-1, seq);
            mpq_class d53 = Ed5(i+1, j-1, seq) + Ed3(i+1, j-1, seq);
            if (seq.V[i][j] + auPenalty(i, j, seq) + constants.multConst[2] + ps.total() <= upper_bound) {
                RNAPartialStructure new_ps(ps);
                new_ps.push(Segment(i, j, lV, seq.V[i][j]));
                new_ps.accumulate(auPenalty(i, j, seq) + constants.multConst[2]);
                pstack.push(new_ps);
                pushed_something = true;
            }
            if (i+1 < j && seq.V[i+1][j] + auPenalty(i+1, j, seq) + constants.multConst[2] + constants.multConst[1] + d5 + ps.total() <= upper_bound) {
                RNAPartialStructure new_ps(ps);
                new_ps.push(Segment(i+1, j, lV, seq.V[i+1][j]));
                new_ps.accumulate(auPenalty(i+1, j, seq) + constants.multConst[2] + constants.multConst[1] + d5);
                new_ps.mark_d5(i);
                pstack.push(new_ps);
                pushed_something = true;
            }
            if (i < j-1 && seq.V[i][j-1] + auPenalty(i, j-1, seq) + constants.multConst[2] + constants.multConst[1] + d3 + ps.total() <= upper_bound) {
                RNAPartialStructure new_ps(ps);
                new_ps.push(Segment(i, j-1, lV, seq.V[i][j-1]));
                new_ps.accumulate(auPenalty(i, j-1, seq) + constants.multConst[2] + constants.multConst[1] + d3);
                new_ps.mark_d3(j);
                pstack.push(new_ps);
                pushed_something = true;
            }
            if (i+1 < j-1 && seq.V[i+1][j-1] + auPenalty(i+1, j-1, seq) + constants.multConst[2] + 2*constants.multConst[1] + d53 + ps.total() <= upper_bound) {
                RNAPartialStructure new_ps(ps);
                new_ps.push(Segment(i+1, j-1, lV, seq.V[i+1][j-1]));
                new_ps.accumulate(auPenalty(i+1, j-1, seq) + constants.multConst[2] + 2*constants.multConst[1] + d53);
                new_ps.mark_d5(i);
                new_ps.mark_d3(j);
                pstack.push(new_ps);
                pushed_something = true;
            }
            break;
        }

        case BOTH_DANGLE:
        {
            mpq_class bonus = Ed5(i, j, seq) + Ed3(i, j, seq) + auPenalty(i, j, seq) + constants.multConst[2];
            if (seq.V[i][j] + bonus + ps.total() <= upper_bound) {
                RNAPartialStructure new_ps(ps);
                new_ps.push(Segment(i, j, lV, seq.V[i][j]));
                new_ps.accumulate(bonus);
                pstack.push(new_ps);
                pushed_something = true;
            }
            break;
        }

        default:
        {
            exit(EXIT_FAILURE);
            break;
        }
        };

        return pushed_something;
    }

    bool NNTM::subopt_traceM(int i, int j, const RNASequenceWithTables& seq, RNAPartialStructure& ps, PartialStructureStack& pstack, mpq_class upper_bound) const {
        // Input specification
        assert (0 <= i);
        assert (i < j);
        assert (j < seq.len());

        bool pushed_something = false;

        if (seq.FM[i][j-1] + constants.multConst[1] + ps.total() <= upper_bound) {
            RNAPartialStructure new_ps(ps);
            new_ps.push(Segment(i, j-1, lM, seq.FM[i][j-1]));
            new_ps.accumulate(constants.multConst[1]);
            pstack.push(new_ps);
            pushed_something = true;
        }

        // case that this whole region is a single branch
        switch (dangles) {
        case NO_DANGLE:
        {
            mpq_class bonus = constants.multConst[2] + auPenalty(i, j, seq);
            if (seq.V[i][j] + bonus + ps.total() <= upper_bound) {
                RNAPartialStructure new_ps(ps);
                new_ps.push(Segment(i, j, lV, seq.V[i][j]));
                new_ps.accumulate(bonus);
                pstack.push(new_ps);
                pushed_something = true;
            }
            break;
        }

        case CHOOSE_DANGLE:
        {
            mpq_class d5 = Ed5(i+1, j, seq);
            mpq_class d3 = Ed3(i, j-1, seq);
            mpq_class d53 = Ed5(i+1, j-1, seq) + Ed3(i+1, j-1, seq);
            if (seq.V[i][j] + constants.multConst[2] + auPenalty(i, j, seq) + ps.total() <= upper_bound) {
                RNAPartialStructure new_ps(ps);
                new_ps.push(Segment(i, j, lV, seq.V[i][j]));
                new_ps.accumulate(constants.multConst[2] + auPenalty(i, j, seq));
                pstack.push(new_ps);
                pushed_something = true;
            }
            if (i+1 < j && seq.V[i+1][j] + constants.multConst[2] + constants.multConst[1] + auPenalty(i+1, j, seq) + d5 + ps.total() <= upper_bound) {
                RNAPartialStructure new_ps(ps);
                new_ps.push(Segment(i+1, j, lV, seq.V[i+1][j]));
                new_ps.accumulate(constants.multConst[2] + constants.multConst[1] + auPenalty(i+1, j, seq) + d5);
                new_ps.mark_d5(i);
                pstack.push(new_ps);
                pushed_something = true;
            }
            if (i < j-1 && seq.V[i][j-1] + constants.multConst[2] + constants.multConst[1] + auPenalty(i, j-1, seq) + d3 + ps.total() <= upper_bound) {
                RNAPartialStructure new_ps(ps);
                new_ps.push(Segment(i, j-1, lV, seq.V[i][j-1]));
                new_ps.accumulate(constants.multConst[2] + constants.multConst[1] + auPenalty(i, j-1, seq) + d3);
                new_ps.mark_d3(j);
                pstack.push(new_ps);
                pushed_something = true;
            }
            if (i+1 < j-1 && seq.V[i+1][j-1] + constants.multConst[2] + 2*constants.multConst[1] + auPenalty(i+1, j-1, seq) + d53 + ps.total() <= upper_bound) {
                RNAPartialStructure new_ps(ps);
                new_ps.push(Segment(i+1, j-1, lV, seq.V[i+1][j-1]));
                new_ps.accumulate(constants.multConst[2] + 2*constants.multConst[1] + auPenalty(i+1, j-1, seq) + d53);
                new_ps.mark_d5(i);
                new_ps.mark_d3(j);
                pstack.push(new_ps);
                pushed_something = true;
            }
            break;
        }

        case BOTH_DANGLE:
        {
            mpq_class bonus = Ed5(i, j, seq) + Ed3(i, j, seq) + auPenalty(i, j, seq) + constants.multConst[2];
            if (seq.V[i][j] + bonus + ps.total() <= upper_bound) {
                RNAPartialStructure new_ps(ps);
                new_ps.push(Segment(i, j, lV, seq.V[i][j]));
                new_ps.accumulate(bonus);
                pstack.push(new_ps);
                pushed_something = true;
            }
            break;
        }

        default:
        {
            exit (EXIT_FAILURE);
            break;
        }
        };

        // case that there are multiple branches
        for (int k = i+TURN+1; k <= j-TURN-1; ++k) {
            switch (dangles) {
            case NO_DANGLE:
            {
                mpq_class bonus = constants.multConst[2] + auPenalty(k+1, j, seq);
                if (seq.FM[i][k] + seq.V[k+1][j] + bonus + ps.total() <= upper_bound) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(i, k, lM, seq.FM[i][k]));
                    new_ps.push(Segment(k+1, j, lV, seq.V[k+1][j]));
                    new_ps.accumulate(bonus);
                    pstack.push(new_ps);
                    pushed_something = true;
                }
                break;
            }

            case CHOOSE_DANGLE:
            {
                mpq_class d5 = Ed5(k+2, j, seq);
                mpq_class d3 = Ed3(k+1, j-1, seq);
                mpq_class d53 = Ed5(k+2, j-1, seq) + Ed3(k+2, j-1, seq);
                if (seq.FM[i][k] + seq.V[k+1][j] + constants.multConst[2] + auPenalty(k+1, j, seq) + ps.total() <= upper_bound) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(i, k, lM, seq.FM[i][k]));
                    new_ps.push(Segment(k+1, j, lV, seq.V[k+1][j]));
                    new_ps.accumulate(constants.multConst[2] + auPenalty(k+1, j, seq));
                    pstack.push(new_ps);
                    pushed_something = true;
                }
                if (k+2 <= j-TURN && seq.FM[i][k] + seq.V[k+2][j] + constants.multConst[2] + constants.multConst[1] + auPenalty(k+2, j, seq) + d5 + ps.total() <= upper_bound) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(i, k, lM, seq.FM[i][k]));
                    new_ps.push(Segment(k+2, j, lV, seq.V[k+2][j]));
                    new_ps.accumulate(constants.multConst[2] + constants.multConst[1] + auPenalty(k+2, j, seq) + d5);
                    new_ps.mark_d5(k+1);
                    pstack.push(new_ps);
                    pushed_something = true;
                }
                if (k+1 <= j-1-TURN && seq.FM[i][k] + seq.V[k+1][j-1] + constants.multConst[2] + constants.multConst[1] + auPenalty(k+1, j-1, seq) + d3 + ps.total() <= upper_bound) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(i, k, lM, seq.FM[i][k]));
                    new_ps.push(Segment(k+1, j-1, lV, seq.V[k+1][j-1]));
                    new_ps.accumulate(constants.multConst[2] + constants.multConst[1] + auPenalty(k+1, j-1, seq) + d3);
                    new_ps.mark_d3(j);
                    pstack.push(new_ps);
                    pushed_something = true;
                }
                if (k+2 <= j-1-TURN && seq.FM[i][k] + seq.V[k+2][j-1] + constants.multConst[2] + 2*constants.multConst[1] + auPenalty(k+2, j-1, seq) + d53 + ps.total() <= upper_bound) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(i, k, lM, seq.FM[i][k]));
                    new_ps.push(Segment(k+2, j-1, lV, seq.V[k+2][j-1]));
                    new_ps.accumulate(constants.multConst[2] + 2*constants.multConst[1] + auPenalty(k+2, j-1, seq) + d53);
                    new_ps.mark_d5(k+1);
                    new_ps.mark_d3(j);
                    pstack.push(new_ps);
                    pushed_something = true;
                }
                break;
            }

            case BOTH_DANGLE:
            {
                mpq_class bonus = Ed5(k+1, j, seq) + Ed3(k+1, j, seq) + constants.multConst[2] + auPenalty(k+1, j, seq);
                if (seq.FM[i][k] + seq.V[k+1][j] + bonus + ps.total() <= upper_bound) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(i, k, lM, seq.FM[i][k]));
                    new_ps.push(Segment(k+1, j, lV, seq.V[k+1][j]));
                    new_ps.accumulate(bonus);
                    pstack.push(new_ps);
                    pushed_something = true;
                }
                break;
            }

            default:
                exit(EXIT_FAILURE);
                break;
            };
        }

        // case that there is a single branch, preceded by free bases ending at position k
        for (int k = i; k <= j-TURN-1; ++k) {

            mpq_class bonus = 0;

            switch (dangles) {
            case NO_DANGLE:
            {
                bonus = constants.multConst[2] + constants.multConst[1]*(k-i+1) + auPenalty(k+1, j, seq);
                if (seq.V[k+1][j] + bonus + ps.total() <= upper_bound) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(k+1, j, lV, seq.V[k+1][j]));
                    new_ps.accumulate(bonus);
                    pstack.push(new_ps);
                    pushed_something = true;
                }
                break;
            }

            case CHOOSE_DANGLE:
            {
                mpq_class d5 = Ed5(k+2, j, seq);
                mpq_class d3 = Ed3(k+1, j-1, seq);
                mpq_class d53 = Ed5(k+2, j-1, seq) + Ed3(k+2, j-1, seq);
                if (seq.V[k+1][j] + constants.multConst[2] + constants.multConst[1]*(k+1 - i) + auPenalty(k+1, j, seq) + ps.total() <= upper_bound) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(k+1, j, lV, seq.V[k+1][j]));
                    new_ps.accumulate(constants.multConst[2] + constants.multConst[1]*(k+1 - i) + auPenalty(k+1, j, seq));
                    pstack.push(new_ps);
                    pushed_something = true;
                }
                if (k+2 <= j-TURN && seq.V[k+2][j] + constants.multConst[2] + constants.multConst[1]*(k+2 - i) + auPenalty(k+2, j, seq) + d5 + ps.total() <= upper_bound) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(k+2, j, lV, seq.V[k+2][j]));
                    new_ps.accumulate(constants.multConst[2] + constants.multConst[1]*(k+2 - i) + auPenalty(k+2, j, seq) + d5);
                    new_ps.mark_d5(k+1);
                    pstack.push(new_ps);
                    pushed_something = true;
                }
                if (k+1 <= j-1-TURN && seq.V[k+1][j-1] + constants.multConst[2] + constants.multConst[1]*(k+1 - i + 1) + auPenalty(k+1, j-1, seq) + d3 + ps.total() <= upper_bound) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(k+1, j-1, lV, seq.V[k+1][j-1]));
                    new_ps.accumulate(constants.multConst[2] + constants.multConst[1]*(k+1 - i + 1) + auPenalty(k+1, j-1, seq) + d3);
                    new_ps.mark_d3(j);
                    pstack.push(new_ps);
                    pushed_something = true;
                }
                if (k+2 <= j-1-TURN && seq.V[k+2][j-1] + constants.multConst[2] + constants.multConst[1]*(k+2 - i + 1) + auPenalty(k+2, j-1, seq) + d53 + ps.total() <= upper_bound) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(k+2, j-1, lV, seq.V[k+2][j-1]));
                    new_ps.accumulate(constants.multConst[2] + constants.multConst[1]*(k+2 - i + 1) + auPenalty(k+2, j-1, seq) + d53);
                    new_ps.mark_d5(k+1);
                    new_ps.mark_d3(j);
                    pstack.push(new_ps);
                    pushed_something = true;
                }
                break;
            }

            case BOTH_DANGLE:
            {
                mpq_class bonus = Ed5(k+1, j, seq) + Ed3(k+1, j, seq) + constants.multConst[2] + constants.multConst[1]*(k-i+1) + auPenalty(k+1, j, seq);
                if (seq.V[k+1][j] + bonus + ps.total() <= upper_bound) {
                    RNAPartialStructure new_ps(ps);
                    new_ps.push(Segment(k+1, j, lV, seq.V[k+1][j]));
                    new_ps.accumulate(bonus);
                    pstack.push(new_ps);
                    pushed_something = true;
                }
                break;
            }

            default:
            {
                exit(EXIT_FAILURE);
                break;
            }
            };
        };

        return pushed_something;
    };
}
