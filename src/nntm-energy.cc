#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "nntm.h"
#include "nndb_constants.h"

#include <gmpxx.h>
#include "boost/multi_array.hpp"

#include <vector>

#ifdef _OPENMP
#include "omp.h"
#endif

namespace pmfe {
    NNTM::NNTM(NNDBConstants constants, dangle_mode dangles):
        constants(constants), dangles(dangles) {
    };

    RNASequenceWithTables NNTM::energy_tables(const RNASequence& inseq) const {
        RNASequenceWithTables seq(inseq, constants.INFINITY_);

        int b, i, j;
        for (b = TURN+1; b <= seq.len() - 1; ++b) {
#ifdef _OPENMP
#pragma omp parallel for private (i,j) schedule(guided)
#endif
            for (i = 1; i <= seq.len() - b; ++i) {
                j = i + b;

                if (seq.can_pair(i, j)) {
                    mpq_class eh = eH(i, j, seq);

                    mpq_class es = eS(i, j, seq) + seq.V[i+1][j-1];

                    seq.VBI[i][j] = calcVBI(i, j, seq);

                    // Multi Loop BEGIN
                    mpq_class d3, d5;
                    d3 = Ed3(i, j, j-1, seq);
                    d5 = Ed5(i, j, i+1, seq);

                    if (dangles == BOTH_DANGLE) { // -d2
                        std::vector<mpq_class> vals;
                        vals.push_back(seq.VM[i][j]);
                        vals.push_back(seq.WMPrime[i+1][j-1] + d3 + d5 + auPenalty(i, j, seq) + constants.multConst[0] + constants.multConst[2]);
                        seq.VM[i][j] = *std::min_element(vals.begin(), vals.end());
                    } else if (dangles == NO_DANGLE) { // -d0
                        std::vector<mpq_class> vals;
                        vals.push_back(seq.VM[i][j]);
                        vals.push_back(seq.WMPrime[i+1][j-1] + auPenalty(i, j, seq) + constants.multConst[0] + constants.multConst[2]);
                        seq.VM[i][j] = *std::min_element(vals.begin(), vals.end());
                    } else { // default
                        std::vector<mpq_class> vals;
                        vals.push_back(seq.VM[i][j]);
                        vals.push_back(seq.WMPrime[i+1][j-1] + auPenalty(i, j, seq) + constants.multConst[0] + constants.multConst[2]);
                        vals.push_back(seq.WMPrime[i+2][j-1] + d5 + auPenalty(i, j, seq) + constants.multConst[0] + constants.multConst[2] + constants.multConst[1]);
                        vals.push_back(seq.WMPrime[i+1][j-2] + d3 + auPenalty(i, j, seq) + constants.multConst[0] + constants.multConst[2] + constants.multConst[1]);
                        vals.push_back(seq.WMPrime[i+2][j-2] + d3 + d5 + auPenalty(i, j, seq) + constants.multConst[0] + constants.multConst[2] + 2*constants.multConst[1]);
                        seq.VM[i][j] = *std::min_element(vals.begin(), vals.end());
                    }

                    // Multi Loop END

                    std::vector<mpq_class> vals;
                    vals.push_back(eh);
                    vals.push_back(es);
                    vals.push_back(seq.VBI[i][j]);
                    vals.push_back(seq.VM[i][j]);
                    seq.V[i][j] = *std::min_element(vals.begin(), vals.end());
                }
                else {
                    seq.V[i][j] = constants.INFINITY_;
                }

                int h;
                for (h = i+TURN+1 ; h <= j-TURN-2; h++) {
                    std::vector<mpq_class> vals;
                    vals.push_back(seq.WMPrime[i][j]);
                    vals.push_back(seq.WM[i][h] + seq.WM[h+1][j]);
                    seq.WMPrime[i][j] = *std::min_element(vals.begin(), vals.end());

                }

                // WM begin
                mpq_class newWM = constants.INFINITY_;

                { // Scoped to protect variable vals
                    std::vector<mpq_class> vals;
                    vals.push_back(newWM);
                    vals.push_back(seq.WMPrime[i][j]);
                    newWM = *std::min_element(vals.begin(), vals.end());
                }

                if (dangles == BOTH_DANGLE) {
                    mpq_class energy = seq.V[i][j] + auPenalty(i, j, seq) + constants.multConst[2];
                    if (i == 1) {
                        energy += Ed3(j, i, seq.len(), seq);
                    } else {
                        energy += Ed3(j, i, i-1, seq);
                    }

                    if (j < seq.len())
                        energy += Ed5(j, i, j+1, seq);

                    std::vector<mpq_class> vals;
                    vals.push_back(energy);
                    vals.push_back(newWM);
                    newWM = *std::min_element(vals.begin(), vals.end());
                } else if (dangles == NO_DANGLE) {
                    std::vector<mpq_class> vals;
                    vals.push_back(seq.V[i][j] + auPenalty(i, j, seq) + constants.multConst[2]);
                    vals.push_back(newWM);
                    newWM = *std::min_element(vals.begin(), vals.end());
                } else { // default
                    std::vector<mpq_class> vals;
                    vals.push_back(newWM);
                    vals.push_back(seq.V[i][j] + auPenalty(i, j, seq) + constants.multConst[2]); // no dangle

                    vals.push_back(seq.V[i+1][j] + Ed3(j, i+1, i, seq) + auPenalty(i+1, j, seq) + constants.multConst[2] + constants.multConst[1]); //i dangle

                    vals.push_back(seq.V[i][j-1] + Ed5(j-1, i, j, seq) + auPenalty(i, j-1, seq) + constants.multConst[2] + constants.multConst[1]);  //j dangle

                    vals.push_back(seq.V[i+1][j-1] + Ed3(j-1, i+1, i, seq) + Ed5(j-1, i+1, j, seq) + auPenalty(i+1, j-1, seq) + constants.multConst[2] + 2*constants.multConst[1]); //i,j dangle

                    newWM = *std::min_element(vals.begin(), vals.end());
                }

                std::vector<mpq_class> vals;
                vals.push_back(newWM);

                vals.push_back(seq.WM[i+1][j] + constants.multConst[1]); //i dangle

                vals.push_back(seq.WM[i][j-1] + constants.multConst[1]); //j dangle

                newWM = *std::min_element(vals.begin(), vals.end());

                seq.WM[i][j] = newWM;
                // WM end
            }
        }

        seq.W[0] = 0;
        for (j = 1; j <= seq.len(); j++) {
            int i;
            mpq_class Wj, Widjd, Wijd, Widj, Wij, Wim1;
            Wj = 0;
            for (i = 1; i < j-TURN; i++) {
                Wij = Widjd = Wijd = Widj = constants.INFINITY_;
                Wim1 = std::min(mpq_class(0), seq.W[i-1]);

                if (dangles == BOTH_DANGLE) { // -d2 option
                    mpq_class energy = seq.V[i][j] + auPenalty(i, j, seq) + Wim1;
                    if (i>1) energy +=  Ed3(j, i, i-1, seq);
                    if (j<seq.len()) energy += Ed5(j, i, j+1, seq);
                    Widjd = energy;

                    std::vector<mpq_class> vals;
                    vals.push_back(Wij);
                    vals.push_back(Widjd);

                    Wij = *std::min_element(vals.begin(), vals.end());
                } else if (dangles == NO_DANGLE) { // -d0 option
                    Wij = seq.V[i][j] + auPenalty(i, j, seq) + Wim1;
                } else { // default
                    Wij = seq.V[i][j] + auPenalty(i, j, seq) + Wim1;
                    Widj = seq.V[i+1][j] + auPenalty(i+1, j, seq) + Ed3(j, i + 1, i, seq) + Wim1;

                    Wijd = seq.V[i][j-1] + auPenalty(i, j-1, seq) + Ed5(j-1, i, j, seq) + Wim1;

                    Widjd = seq.V[i+1][j-1] + auPenalty(i+1, j-1, seq) + Ed3(j-1, i+1, i, seq) + Ed5(j-1, i+1, j, seq) + Wim1;

                    std::vector<mpq_class> vals;
                    vals.push_back(Wij);
                    vals.push_back(Widj);
                    vals.push_back(Wijd);
                    vals.push_back(Widjd);

                    Wij = *std::min_element(vals.begin(), vals.end());
                }

                std::vector<mpq_class> vals;
                vals.push_back(Wj);
                vals.push_back(Wij);

                Wj = *std::min_element(vals.begin(), vals.end());
            }

            std::vector<mpq_class> vals;
            vals.push_back(Wj);
            vals.push_back(seq.W[j-1]);

            seq.W[j] = *std::min_element(vals.begin(), vals.end());
        }

        return seq;
    };

    mpq_class NNTM::minimum_energy(const RNASequenceWithTables& seq) const {
        return seq.W[seq.len()];
    };


    // TODO: Replace these with semantic Ed3 and Ed5
    mpq_class NNTM::Ed3(int i, int j, int k, const RNASequenceWithTables& seq) const {
        return constants.dangle[seq.base(i)][seq.base(j)][seq.base(k)][1];
    }

    mpq_class NNTM::Ed5(int i, int j, int k, const RNASequenceWithTables& seq) const {
        return constants.dangle[seq.base(i)][seq.base(j)][seq.base(k)][0];
    }

    mpq_class NNTM::auPenalty(int i, int j, const RNASequenceWithTables& seq) const {
        // Return the penalty for non-GC pairs
        int base_i = seq.base(i);
        int base_j = seq.base(j);
        if (
            (base_i == BASE_U and (base_j == BASE_A or base_j == BASE_G)) or
            (base_j == BASE_U and (base_i == BASE_A or base_i == BASE_G))
            ) {
            return constants.auend;
        } else {
            return 0;
        }
    }

    mpq_class NNTM::eL(int i, int j, int ip, int jp, const RNASequenceWithTables& seq) const {
        mpq_class energy = constants.INFINITY_;
        mpq_class loginc = 0;

        /*SH: These calculations used to incorrectly be within the bulge loop code, moved out here. */
        int size1 = ip - i - 1;
        int size2 = j - jp - 1;
        int size = size1 + size2;

        int lopsided = abs(size1 - size2); /* define the asymmetry of an interior loop */

        int pindex = std::min(2, std::min(size1, size2));
        mpq_class lvalue = lopsided * constants.poppen[pindex];
        mpq_class minterm = std::min(constants.maxpen, lvalue);

        if (size1 == 0 || size2 == 0) {
            if (size > 30) {
                /* AM: Does not depend upon i and j and ip and jp - Stacking Energies */
                loginc = constants.prelog * log((double) size / 30.0); // Taking a log and then casting to rational--gross! But we gots ta do what we gots ta do
                energy = constants.bulge[30] + constants.eparam[2] + loginc + auPenalty(i, j, seq)
                    + auPenalty(ip, jp, seq);
            } else if (size <= 30 && size != 1) {
                /* Does not depend upon i and j and ip and jp - Stacking Energies  */
                energy = constants.bulge[size] + constants.eparam[2];
                energy += auPenalty(i, j, seq) + auPenalty(ip, jp, seq);
            } else if (size == 1) {
                energy = constants.stack[seq.base(i)][seq.base(j)][seq.base(ip)][seq.base(jp)]
                    + constants.bulge[size] + constants.eparam[2];
            }
        } else {
            /* Internal loop */

            if (size > 30) {
                loginc = constants.prelog * log((double) size / 30.0);

                if (!((size1 == 1 || size2 == 1) && constants.gail)) { /* normal internal loop with size > 30*/
                    energy = constants.tstki[seq.base(i)][seq.base(j)][seq.base(i + 1)][seq.base(j - 1)] +
                        constants.tstki[seq.base(jp)][seq.base(ip)][seq.base(jp + 1)][seq.base(ip - 1)] + constants.inter[30] + loginc +
                        constants.eparam[3] + minterm;
                } else { /* if size is more than 30 and it is a grossely asymmetric internal loop and gail is not 0*/
                    energy = constants.tstki[seq.base(i)][seq.base(j)][BASE_A][BASE_A] + constants.tstki[seq.base(jp)][seq.base(ip)][BASE_A][BASE_A]
                        + constants.inter[30] + loginc + constants.eparam[3] + minterm;
                }
            }
            else if (size1 == 2 && size2 == 2) { /* 2x2 internal loop */
                energy = constants.iloop22[seq.base(i)][seq.base(ip)][seq.base(j)][seq.base(jp)][seq.base(i+1)][seq.base(i+2)][seq.base(j-1)][seq.base(j-2)];
            } else if (size1 == 1 && size2 == 2) {
                energy = constants.iloop21[seq.base(i)][seq.base(j)][seq.base(i + 1)][seq.base(j - 1)][seq.base(j - 2)][seq.base(ip)][seq.base(jp)];
            } else if (size1 == 2 && size2 == 1) { /* 1x2 internal loop */
                energy = constants.iloop21[seq.base(jp)][seq.base(ip)][seq.base(j - 1)][seq.base(i + 2)][seq.base(i + 1)][seq.base(j)][seq.base(i)];
            } else if (size == 2) { /* 1*1 internal loops */
                energy = constants.iloop11[seq.base(i)][seq.base(i + 1)][seq.base(ip)][seq.base(j)][seq.base(j - 1)][seq.base(jp)];
            }
            //else if ((size1 == 2 && size2 == 3) || (size1 == 3 && size2 == 2)) {
            //	return constants.tstacki23[seq.base(i)][seq.base(j)][seq.base(i + 1)][seq.base(j - 1)] +
            //		constants.tstacki23[seq.base(jp)][seq.base(ip)][seq.base(jp + 1)][seq.base(ip - 1)];
            //}
            else if ((size1 == 1 || size2 == 1) && constants.gail) { /* gail = (Grossly Asymmetric Interior Loop Rule) (on/off <-> 1/0)  */
                energy = constants.tstki[seq.base(i)][seq.base(j)][BASE_A][BASE_A] + constants.tstki[seq.base(jp)][seq.base(ip)][BASE_A][BASE_A]
                    + constants.inter[size] + loginc + constants.eparam[3] + minterm;
            } else { /* General Internal loops */
                energy = constants.tstki[seq.base(i)][seq.base(j)][seq.base(i + 1)][seq.base(j - 1)] + constants.tstki[seq.base(jp)][seq.base(ip)][seq.base(jp + 1)][seq.base(ip - 1)] + constants.inter[size]
                    + loginc + constants.eparam[3] + minterm;
            }
        }

        return energy;
    }

    mpq_class NNTM::eH(int i, int j, const RNASequenceWithTables& seq) const {
        /*  Hairpin loop for all the bases between i and j */
        /*  size for size of the loop, energy is the result, loginc is for the extrapolation for loops bigger than 30 */
        int size;
        mpq_class loginc;
        mpq_class energy;
        mpq_class tlink;

        energy = constants.INFINITY_;

        size = j - i - 1; /*  size is the number of bases in the loop, when the closing pair is excluded */

        /*  look in hairpin, and be careful that there is only 30 values */

        if (size > 30) {
            loginc = constants.prelog * log(((double) size) / 30.0);
            energy = constants.hairpin[30] + loginc + constants.tstkh[seq.base(i)][seq.base(j)][seq.base(i + 1)][seq.base(j - 1)] + constants.eparam[4]; /* size penalty + terminal mismatch stacking energy*/
        }

        else if (size <= 30 && size > 4) {
            energy = constants.hairpin[size] + constants.tstkh[seq.base(i)][seq.base(j)][seq.base(i + 1)][seq.base(j - 1)] + constants.eparam[4]; /* size penalty + terminal mismatch stacking energy*/
        }

        else if (size == 4) {
            std::string loopkey = seq.subsequence(i, i+6);
            mpq_class tlink = 0; // Loop contribution is typically 0
            if (constants.tloop.count(loopkey) != 0) {
                tlink = constants.tloop.find(loopkey)->second; // But some loops have special contributions, stored in this table
            }
            energy = tlink + constants.hairpin[size] + constants.tstkh[seq.base(i)][seq.base(j)][seq.base(i + 1)][seq.base(j - 1)] + constants.eparam[4];
        }

        else if (size == 3) {
            /*  triloop... For the moment, the file triloop.dat is empty */
            /*  else, should have a treatment like the one if size==4 */
            energy = constants.hairpin[size];
            /* AM: Don't include stacking energy terms for triloopls */
            /* + constants.tstkh[seq.base(i)][seq.base(j)][seq.base(i+1)][seq.base(j-1)]  */
            /* + constants.eparam[4]; */
            /*  Must be another penalty for terminal AU... Not sure of this */
            energy += auPenalty(i, j, seq);
        }

        else if (size < 3 && size != 0) {
            /*  no terminal mismatch */
            energy = constants.hairpin[size] + constants.eparam[4];
        } else if (size == 0)
            energy = constants.INFINITY_;

        /*  GGG Bonus => GU closure preceded by GG */
        /*  i-2 = i-1 = i = G, and j = U; i < j */
        if (i > 2) {
            if (seq.base(i - 2) == BASE_G && seq.base(i - 1) == BASE_G && seq.base(i) == BASE_G
                && seq.base(j) == BASE_U) {
                energy += constants.gubonus;
                /*  printf ("\n GGG bonus for i %d j %d ", i, j); */
            }
        }

        /*  Poly-C loop => How many C are needed for being a poly-C loop */
        tlink = 1;
        for (int index = 1; (index <= size) && (tlink == 1); ++index) {
            if (seq.base(i + index) != BASE_C)
                tlink = 0;
        }

        if (tlink == 1) {
            if (size == 3) {
                energy += constants.c3;
            } else {
                energy += constants.cint + size * constants.cslope;
            }
        }

        return energy;
    }

    mpq_class NNTM::eS(int i, int j, const RNASequenceWithTables& seq) const {
        /*  not sure about constants.eparam[1], come from MFold.. = 0 */
        return constants.stack[seq.base(i)][seq.base(j)][seq.base(i+1)][seq.base(j-1)] + constants.eparam[1];
    }

    mpq_class NNTM::calcVBI(int i, int j, const RNASequenceWithTables& seq) const {
        mpq_class VBIij = constants.INFINITY_;

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

            for (int q = minq; q <= maxq; q++) {
                if (q - p > TURN && seq.can_pair(p, q)) {
                    std::vector<mpq_class> vals;
                    vals.push_back(eL(i, j, p, q, seq) + seq.V[p][q]);
                    vals.push_back(VBIij);
                    VBIij = *std::min_element(vals.begin(), vals.end());
                }
            }
        }

        return VBIij;
    }
};
