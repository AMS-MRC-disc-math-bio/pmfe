// Copyright (c) 2015 Andrew Gainer-Dewar.

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdexcept>

#include "nntm.h"
#include "nndb_constants.h"

#include <gmpxx.h>
#include <boost/multi_array.hpp>

#include <vector>

namespace pmfe {
    RNAStructureWithScore NNTM::mfe_structure(const RNASequenceWithTables& seq) const {
        RNAStructure structure(seq);
        ScoreVector score;

        traceW(seq.len()-1, seq, structure, score);

        ScoreVector newscore = this->score(structure);

        if (newscore.energy != score.energy) {
            throw std::logic_error("Energy calculation was inconsistent!");
        }

        RNAStructureWithScore result (structure, newscore);
        return result;
    }

    void NNTM::traceW(int j, const RNASequenceWithTables& seq, RNAStructure& structure, ScoreVector& score) const {
        bool finished = false;
        mpq_class wim1;

        if (j <= 0)
            return;

        for (int i = 0; i < j && !finished; i++) {
            if (j-i < TURN) continue;

            if (i > 0) {
                wim1 = std::min(mpq_class(0), seq.W[i-1]);
            } else {
                wim1 = 0;
            }

            if (dangles == BOTH_DANGLE) {
                mpq_class e_dangles = 0;
                if (i > 0) {
                    e_dangles += Ed5(i, j, seq);
                }

                if (j < seq.len() - 1) {
                    e_dangles += Ed3(i, j, seq);
                }

                if (seq.W[j] == seq.V[i][j] + auPenalty(i, j, seq) + e_dangles + wim1) {
                    finished = true;
                    score.energy += (auPenalty(i, j, seq) + e_dangles);
                    traceV(i, j, seq, structure, score);
                    traceW(i-1, seq, structure, score);
                    break;
                };
            } else if (dangles == NO_DANGLE) {
                if (seq.W[j] == seq.V[i][j] + auPenalty(i, j, seq) + wim1) {
                    finished = true;
                    score.energy += auPenalty(i, j, seq);
                    traceV(i, j, seq, structure, score);
                    traceW(i-1, seq, structure, score);
                    break;
                };
            } else { // default
                if (seq.W[j] == seq.V[i][j] + auPenalty(i, j, seq) + wim1) {
                    finished = true;
                    score.energy += auPenalty(i, j, seq);
                    traceV(i, j, seq, structure, score);
                    traceW(i-1, seq, structure, score);
                    break;
                } else if (seq.W[j] ==  seq.V[i][j-1] + auPenalty(i, j-1, seq) + Ed3(i, j-1, seq) + wim1) {
                    finished = true;
                    score.energy += (auPenalty(i, j-1, seq) + Ed3(i, j-1, seq));
                    structure.mark_d3(j);
                    traceV(i, j-1, seq, structure, score);
                    traceW(i-1, seq, structure, score);
                    break;
                } else if (seq.W[j] == seq.V[i+1][j] + auPenalty(i+1, j, seq) + Ed5(i+1, j, seq) + wim1){
                    finished = true;
                    score.energy += (auPenalty(i+1, j, seq) + Ed5(i+1, j, seq));
                    structure.mark_d5(i);
                    traceV(i + 1, j, seq, structure, score);
                    traceW(i-1, seq, structure, score);
                    break;
                } else if (seq.W[j] == seq.V[i+1][j-1] + auPenalty(i+1, j-1, seq) + Ed5(i+1, j-1, seq) + Ed3(i+1, j-1, seq) + wim1) {
                    finished = true;
                    score.energy += (auPenalty(i+1, j-1, seq) + Ed5(i+1, j-1, seq) + Ed3(i+1, j-1, seq));
                    structure.mark_d3(j);
                    structure.mark_d5(i);
                    traceV(i+1, j-1, seq, structure, score);
                    traceW(i-1, seq, structure, score);
                    break;
                }
            }
        }

        if (seq.W[j] == seq.W[j-1] && !finished) traceW(j-1, seq, structure, score);
        return;
    }

    mpq_class NNTM::traceV(int i, int j, const RNASequenceWithTables& seq, RNAStructure& structure, ScoreVector& score) const {
        mpq_class a, b, c, d;
        mpq_class Vij;
        if (j-i < TURN)  return constants.INFINITY_;

        // TODO: Eliminate silly intermediate variables
        a = eH(i, j, seq);

        b = eS(i, j, seq) + seq.V[i + 1][j - 1];
        c = seq.VBI[i][j];
        d = seq.VM[i][j];

        Vij = seq.V[i][j];
        structure.mark_pair(i, j);

        if (Vij == a ) {
            score.energy += eH(i, j, seq);
            return Vij;
        } else if (Vij == b) {
            score.energy += eS(i, j, seq);
            traceV(i+1, j-1, seq, structure, score);
            return Vij;
        } else if (Vij == c) {
            traceVBI(i, j, seq, structure, score);
            return Vij;
        } else if (Vij == d) {
            mpq_class eVM = traceVM(i, j, seq, structure, score);
            score.energy += (Vij-eVM);
            return Vij;
        }

        return 0;
    }

    mpq_class NNTM::traceVBI(int i, int j, const RNASequenceWithTables& seq, RNAStructure& structure, ScoreVector& score) const {
        mpq_class VBIij;
        int ip, jp;
        int ifinal, jfinal;

        ifinal = 0;
        jfinal = 0;

        for (ip = i + 1; ip < j - 1; ip++) {
            for (jp = ip + 1; jp < j; jp++) {
                VBIij = eL(i, j, ip, jp, seq) + seq.V[ip][jp];
                if (VBIij == seq.VBI[i][j]){
                    ifinal = ip;
                    jfinal = jp;
                    break;
                }
            }
            if (jp != j) break;
        }

        score.energy += eL(i, j, ifinal, jfinal, seq);

        return traceV(ifinal, jfinal, seq, structure, score);
    }

    mpq_class NNTM::traceVM(int i, int j, const RNASequenceWithTables& seq, RNAStructure& structure, ScoreVector& score) const {
        mpq_class eVM = 0;

        if (dangles == BOTH_DANGLE) {
            if (seq.VM[i][j] == seq.WMPrime[i+1][j-1] + constants.multConst[0] + constants.multConst[2] + auPenalty(i, j, seq) + Ed5(i, j, seq, true) + Ed3(i, j, seq, true)) {
                eVM += traceWMPrime(i+1, j-1, seq, structure, score);
                score.multiloops++;
                score.branches++;
            }
        } else if (dangles == NO_DANGLE) {
            if (seq.VM[i][j] == seq.WMPrime[i+1][j-1] + constants.multConst[0] + constants.multConst[2] + auPenalty(i, j, seq) ) {
                eVM += traceWMPrime(i+1, j-1, seq, structure, score);
                score.multiloops++;
                score.branches++;
            }
        } else {
            if (seq.VM[i][j] == seq.WMPrime[i+1][j-1] + constants.multConst[0] + constants.multConst[2] + auPenalty(i, j, seq) ) {
                eVM += traceWMPrime(i+1, j-1, seq, structure, score);
                score.multiloops++;
                score.branches++;
            } else if (seq.VM[i][j] == seq.WMPrime[i+2][j-1] + constants.multConst[0] + constants.multConst[2] + auPenalty(i, j, seq) + Ed5(i, j, seq, true) + constants.multConst[1]) {
                eVM += traceWMPrime(i+2, j-1, seq, structure, score);
                structure.mark_d3(i+1);
                score.multiloops++;
                score.branches++;
                score.unpaired++;
            }
            else if (seq.VM[i][j] == seq.WMPrime[i+1][j-2] + constants.multConst[0] + constants.multConst[2] + auPenalty(i, j, seq) + Ed3(i, j, seq, true) + constants.multConst[1]) {
                eVM += traceWMPrime(i+1, j-2, seq, structure, score);
                structure.mark_d5(j-1);
                score.multiloops++;
                score.branches++;
                score.unpaired++;
            } else if (seq.V[i][j] ==  seq.WMPrime[i+2][j-2] + constants.multConst[0] + constants.multConst[2] + auPenalty(i, j, seq) + Ed5(i, j, seq, true) + Ed3(i, j, seq, true) + constants.multConst[1]*2) {
                eVM += traceWMPrime(i+2, j-2, seq, structure, score);
                structure.mark_d3(i+1);
                structure.mark_d5(j-1);
                score.multiloops++;
                score.branches++;
                score.unpaired += 2;
            }
        }

        return eVM;
    }

    mpq_class NNTM::traceWMPrime(int i, int j, const RNASequenceWithTables& seq, RNAStructure& structure, ScoreVector& score) const {
        int done=0, h;
        mpq_class energy = 0;

        for (h = i; h < j && !done; h++) {
            if (seq.WM[i][h] + seq.WM[h+1][j] == seq.WMPrime[i][j]) {
                energy += traceWM(i, h, seq, structure, score);
                energy += traceWM(h+1, j, seq, structure, score);
                done = 1;
                break;
            }
        }
        return energy;
    }

    mpq_class NNTM::traceWM(int i, int j, const RNASequenceWithTables& seq, RNAStructure& structure, ScoreVector& score) const {
        assert(i < j);
        int done = 0;
        mpq_class eWM = 0;

        if (!done && seq.WM[i][j] == seq.WMPrime[i][j]) {
            eWM += traceWMPrime(i, j, seq, structure, score);
            done = 1;
        }

        if (!done){
            if (dangles == BOTH_DANGLE) {
                if (seq.WM[i][j] == seq.V[i][j] + auPenalty(i, j, seq) + constants.multConst[2] + Ed5(i, j, seq) + Ed3(i, j, seq)) {
                    eWM += traceV(i, j, seq, structure, score);
                    score.branches++;
                    done = 1;
                }
            } else if (dangles == NO_DANGLE) {
                if (seq.WM[i][j] == seq.V[i][j] + auPenalty(i, j, seq) + constants.multConst[2]) {
                    eWM += traceV(i, j, seq, structure, score);
                    score.branches++;
                    done = 1;
                }
            } else  {
                if (seq.WM[i][j] == seq.V[i][j] + auPenalty(i, j, seq) + constants.multConst[2]) {
                    eWM += traceV(i, j, seq, structure, score);
                    score.branches++;
                    done = 1;
                } else if (seq.WM[i][j] == seq.V[i+1][j] + Ed5(i+1, j, seq) + auPenalty(i+1, j, seq) + constants.multConst[2] + constants.multConst[1]) {
                    eWM += traceV(i+1, j, seq, structure, score);
                    structure.mark_d5(i);
                    score.branches++;
                    score.unpaired++;
                    done = 1;
                } else if (seq.WM[i][j] == seq.V[i][j-1] + Ed3(i, j-1, seq) + auPenalty(i, j-1, seq) + constants.multConst[2] + constants.multConst[1]) {
                    eWM += traceV(i, j-1, seq, structure, score);
                    structure.mark_d3(j);
                    score.branches++;
                    score.unpaired++;
                    done = 1;
                } else if (seq.WM[i][j] == seq.V[i+1][j-1] + Ed5(i+1, j-1, seq) + Ed3(i+1, j-1, seq) + auPenalty(i+1, j-1, seq) + constants.multConst[2] + constants.multConst[1]*2) {
                    eWM += traceV(i+1, j-1, seq, structure, score);
                    structure.mark_d5(i);
                    structure.mark_d3(j);
                    score.branches++;
                    score.unpaired += 2;
                    done = 1;
                }
            }
        }

        if (!done){
            if (seq.WM[i][j] == seq.WM[i+1][j] + constants.multConst[1]) {
                done = 1;
                eWM += traceWM(i+1, j, seq, structure, score);
                score.unpaired++;
            } else if (seq.WM[i][j] == seq.WM[i][j-1] + constants.multConst[1]) {
                done = 1;
                eWM += traceWM(i, j-1, seq, structure, score);
                score.unpaired++;
            }
        }

        return eWM;
    }
}
