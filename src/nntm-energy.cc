// Copyright (c) 2015 Andrew Gainer-Dewar

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>

#include "nntm.h"
#include "nndb_constants.h"
#include "rational.h"
#include "minbox.h"

#include <boost/bind.hpp>

#define BOOST_LOG_DYN_LINK 1 // Fix an issue with dynamic library loading
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

namespace pmfe {
    NNTM::NNTM(const NNDBConstants& constants, dangle_mode dangles):
        constants(constants),
        dangles(dangles)
    {};

    RNASequenceWithTables NNTM::energy_tables(const RNASequence& inseq) const {
        RNASequenceWithTables seq(inseq);
        populate_energy_tables(seq);
        return seq;
    };

    void NNTM::populate_energy_tables(RNASequenceWithTables& seq) const {
        /*
          Construct the energy tables for the DP algorithm
        */
        // Input specification
        assert(not seq.energy_tables_populated);

        // Populate V, VM, VBI, WM, and WMPrime
        for (int b = TURN+1; b <= seq.len() - 1; ++b) {
#pragma omp parallel for shared(seq)
            for (int i = 0; i <= seq.len() - 1 - b; ++i) {
                populate_energy_tables(i, i+b, seq);
            }
        }

        // Populate W
        for (int j = 0; j <= seq.len() - 1; ++j) {
            if (j <= TURN) {
                seq.W[j] = 0;
                continue;
            }

            MinBox<Rational> w_vals;
            w_vals.insert(Rational::infinity());

            for (int i = 0; i < j-TURN; i++) {
                Rational Wim1;
                if (i > 0) {
                    Wim1 = seq.W[i-1];
                } else {
                    Wim1 = 0;
                }

                switch (dangles) {
                case BOTH_DANGLE:
                    {
                        Rational Widjd = seq.V[i][j] + auPenalty(i, j, seq) + Wim1;
                        if (i > 0) {
                            Widjd += Ed5(i, j, seq);
                        }

                        if (j < seq.len() - 1) {
                            Widjd += Ed3(i, j, seq);
                        }

                        w_vals.insert(Widjd);
                        break;
                    }

                case NO_DANGLE:
                    {
                        w_vals.insert(seq.V[i][j] + auPenalty(i, j, seq) + Wim1);
                        break;
                    }

                case CHOOSE_DANGLE:
                    {
                        w_vals.insert(seq.V[i][j] + auPenalty(i, j, seq) + Wim1);
                        w_vals.insert(seq.V[i+1][j] + auPenalty(i+1, j, seq) + Ed5(i+1, j, seq) + Wim1);
                        w_vals.insert(seq.V[i][j-1] + auPenalty(i, j-1, seq) + Ed3(i, j-1, seq) + Wim1);
                        w_vals.insert(seq.V[i+1][j-1] + auPenalty(i+1, j-1, seq) + Ed5(i+1, j-1, seq) + Ed3(i+1, j-1, seq) + Wim1);
                        break;
                    }

                default:
                    throw std::logic_error("Invalid dangle mode.");
                    break;
                }
            }

            w_vals.insert(seq.W[j-1]); // Base j is free
            w_vals.insert(0); // All bases up to j are free

            seq.W[j] = w_vals.minimum();
        }

        seq.energy_tables_populated = true;
    }

    void NNTM::populate_energy_tables(int i, int j, RNASequenceWithTables& seq) const {
        // Input specification
        assert (0 <= i);
        assert (j < seq.len());
        assert (i < j);

        if (seq.can_pair(i, j)) {
            MinBox<Rational> vm_vals;
            vm_vals.insert(Rational::infinity());

            Rational d3, d5;
            d3 = Ed3(i, j, seq, true);
            d5 = Ed5(i, j, seq, true);

            switch (dangles) {
            case BOTH_DANGLE:
                {
                    vm_vals.insert(seq.WMPrime[i+1][j-1] + d3 + d5 + auPenalty(i, j, seq) + constants.multConst[0] + constants.multConst[2]);
                    break;
                }

            case NO_DANGLE:
                {
                    vm_vals.insert(seq.WMPrime[i+1][j-1] + auPenalty(i, j, seq) + constants.multConst[0] + constants.multConst[2]);
                    break;
                }

            case CHOOSE_DANGLE:
                {
                    vm_vals.insert(seq.WMPrime[i+1][j-1] + auPenalty(i, j, seq) + constants.multConst[0] + constants.multConst[2]);
                    vm_vals.insert(seq.WMPrime[i+2][j-1] + d5 + auPenalty(i, j, seq) + constants.multConst[0] + constants.multConst[2] + constants.multConst[1]);
                    vm_vals.insert(seq.WMPrime[i+1][j-2] + d3 + auPenalty(i, j, seq) + constants.multConst[0] + constants.multConst[2] + constants.multConst[1]);
                    vm_vals.insert(seq.WMPrime[i+2][j-2] + d3 + d5 + auPenalty(i, j, seq) + constants.multConst[0] + constants.multConst[2] + 2*constants.multConst[1]);
                    break;
                }

            default:
                throw std::logic_error("Invalid dangle mode.");
                break;
            }

            seq.VM[i][j] = vm_vals.minimum();

            MinBox<Rational> v_vals;
            v_vals.insert(Rational::infinity());
            v_vals.insert(seq.VM[i][j]);

            v_vals.insert(eH(i, j, seq));
            v_vals.insert(eS(i, j, seq) + seq.V[i+1][j-1]);

            seq.VBI[i][j] = calcVBI(i, j, seq);
            v_vals.insert(seq.VBI[i][j]);

            seq.V[i][j] = v_vals.minimum();
        } else {
            seq.V[i][j] = Rational::infinity();
        }

        MinBox<Rational> wmp_vals;
        wmp_vals.insert(Rational::infinity());

        for (int h = i+TURN+1 ; h <= j-TURN-2; ++h) {
            wmp_vals.insert(seq.WM[i][h] + seq.WM[h+1][j]);
        }

        seq.WMPrime[i][j] = wmp_vals.minimum();

        // WM begin
        MinBox<Rational> wm_vals;
        wm_vals.insert(Rational::infinity());
        wm_vals.insert(seq.WMPrime[i][j]);

        switch (dangles) {
        case BOTH_DANGLE:
            {
                Rational energy = seq.V[i][j] + auPenalty(i, j, seq) + constants.multConst[2];

                if (i > 0) {
                    energy += Ed5(i, j, seq);
                }

                if (j < seq.len() - 1) {
                    energy += Ed3(i, j, seq);
                }

                wm_vals.insert(energy);
                break;
            }

        case NO_DANGLE: {
            wm_vals.insert(seq.V[i][j] + auPenalty(i, j, seq) + constants.multConst[2]);
            break;
        }

        case CHOOSE_DANGLE:
            {
                wm_vals.insert(seq.V[i][j] + auPenalty(i, j, seq) + constants.multConst[2]); // no dangle

                wm_vals.insert(seq.V[i+1][j] + Ed5(i+1, j, seq) + auPenalty(i+1, j, seq) + constants.multConst[2] + constants.multConst[1]); //i dangle

                wm_vals.insert(seq.V[i][j-1] + Ed3(i, j-1, seq) + auPenalty(i, j-1, seq) + constants.multConst[2] + constants.multConst[1]);  //j dangle

                wm_vals.insert(seq.V[i+1][j-1] + Ed5(i+1, j-1, seq) + Ed3(i+1, j-1, seq) + auPenalty(i+1, j-1, seq) + constants.multConst[2] + 2*constants.multConst[1]); //i,j dangle

                break;
            }

        default:
            throw std::logic_error("Invalid dangle mode.");
            break;
        }

        wm_vals.insert(seq.WM[i+1][j] + constants.multConst[1]); //i dangle

        wm_vals.insert(seq.WM[i][j-1] + constants.multConst[1]); //j dangle

        seq.WM[i][j] = wm_vals.minimum();
        // WM end
    }

    Rational NNTM::minimum_energy(RNASequenceWithTables& seq) const {
        /*
          Return the minimum energy of a structure on this sequence
        */
        if (not seq.energy_tables_populated) {
            populate_energy_tables(seq);
        }

        return seq.W[seq.len()-1];
    };


    // dangle on the 5' end of (i, j)
    // if inside==true, dangle i+1 instead of i-1
    Rational NNTM::Ed5(int i, int j, const RNASequence& seq, bool inside) const {
        // Input specification
        assert (i >= 0 and i < seq.len());
        assert (j >= 0 and j < seq.len());

        Rational penalty = 0;

        if (i > 0) {
            if (inside)
                penalty = constants.dangle[seq.base(i)][seq.base(j)][seq.base(i+1)][0];
            else
                penalty = constants.dangle[seq.base(j)][seq.base(i)][seq.base(i-1)][1];
        } else {
            if (inside)
                penalty = constants.dangle[seq.base(i)][seq.base(j)][seq.base(i+1)][0];
            else
                penalty = constants.dangle[seq.base(j)][seq.base(i)][seq.base(seq.len()-1)][1];
        }

        return penalty;
    }

    // dangle on the 3' end of (i, j)
    // if inside==true, dangle j-1 instead of j+1
    Rational NNTM::Ed3(int i, int j, const RNASequence& seq, bool inside) const {
        // Input specification
        assert (i >= 0 and i < seq.len());
        assert (j >= 0 and j < seq.len());

        Rational penalty = 0;

        if (j < seq.len()-1) {
            if (inside)
                penalty = constants.dangle[seq.base(i)][seq.base(j)][seq.base(j-1)][1];
            else
                penalty = constants.dangle[seq.base(j)][seq.base(i)][seq.base(j+1)][0];
        } else {
            if (inside)
                penalty = constants.dangle[seq.base(i)][seq.base(j)][seq.base(j-1)][1];
            else
                penalty = constants.dangle[seq.base(j)][seq.base(i)][seq.base(0)][0];
        }

        return penalty;
    }

    Rational NNTM::auPenalty(int i, int j, const RNASequence& seq) const {
        /*
          Return the pairing penalty for (i, j); this is nonzero unless it is a GC pair
        */
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

    Rational NNTM::eLL(int size) const {
        /*
          Compute the energy correction for a long loop
        */

        // To match the GTMFE algorithm, we force this to have two decimal digits of precision
        // TODO: Clean this up by interfacing with fancified NNDB
        if (constants.params.dummy_scaling == 0) {
            return 0;
        } else {
            Rational prelog = constants.prelog / constants.params.dummy_scaling;
            Rational result = (Integer) (100 * prelog.get_d() * log((double) size / 30.0));
            result *= Rational(constants.params.dummy_scaling / 100);
            return result;
        }
    }

    Rational NNTM::eL(int i, int j, int ip, int jp, const RNASequence& seq) const {
        /*
          Compute the energy of an internal loop between pairs (i, j) and (ip, jp)
        */
        // Input specification
        assert (i >= 0 and i < seq.len());
        assert (j >= 0 and j < seq.len());
        assert (i < ip and ip < jp and jp < j);

        Rational energy = Rational::infinity();

        /*SH: These calculations used to incorrectly be within the bulge loop code, moved out here. */
        int size1 = ip - i - 1;
        int size2 = j - jp - 1;
        int size = size1 + size2;

        int lopsided = abs(size1 - size2); /* define the asymmetry of an interior loop */

        int pindex = std::min(2, std::min(size1, size2));
        Rational lvalue = lopsided * constants.poppen[pindex];
        Rational penterm;

        if (constants.params.dummy_scaling >= 0) {
            penterm = std::min(constants.maxpen, lvalue);
        } else {
            penterm = std::max(constants.maxpen, lvalue);
        }

        if (size1 == 0 or size2 == 0) {
            if (size > 30) {
                /* AM: Does not depend upon i and j and ip and jp - Stacking Energies */
                energy = constants.bulge[30]
                    + eLL(size)
                    + auPenalty(i, j, seq)
                    + auPenalty(ip, jp, seq);
            } else if (size <= 30 and size > 1) {
                /* Does not depend upon i and j and ip and jp - Stacking Energies  */
                energy = constants.bulge[size]
                    + auPenalty(i, j, seq)
                    + auPenalty(ip, jp, seq);
            } else if (size == 1) {
                energy = constants.stack[seq.base(i)][seq.base(j)][seq.base(ip)][seq.base(jp)]
                    + constants.bulge[size];
            } else if (size == 0) {
                // Degenerate case that this is actually a stack, included for ease of coding elsewhere
                energy = constants.stack[seq.base(i)][seq.base(j)][seq.base(ip)][seq.base(jp)];
            }
        } else {
            /* Internal loop */

            if (size > 30) {
                if (not ((size1 == 1 or size2 == 1) and constants.gail)) { /* normal internal loop with size > 30*/
                    energy = constants.tstki[seq.base(i)][seq.base(j)][seq.base(i + 1)][seq.base(j - 1)]
                        + constants.tstki[seq.base(jp)][seq.base(ip)][seq.base(jp + 1)][seq.base(ip - 1)] + constants.inter[30]
                        + eLL(size)
                        + penterm;
                } else { /* if size is more than 30 and it is a grossely asymmetric internal loop and gail is not 0*/
                    energy = constants.tstki[seq.base(i)][seq.base(j)][BASE_A][BASE_A]
                        + constants.tstki[seq.base(jp)][seq.base(ip)][BASE_A][BASE_A]
                        + constants.inter[30]
                        + eLL(size)
                        + penterm;
                }
            } else if (size1 == 2 and size2 == 2) { /* 2x2 internal loop */
                energy = constants.iloop22[seq.base(i)][seq.base(ip)][seq.base(j)][seq.base(jp)][seq.base(i+1)][seq.base(j-1)][seq.base(i+2)][seq.base(j-2)];
            } else if (size1 == 1 and size2 == 2) {
                energy = constants.iloop21[seq.base(i)][seq.base(j)][seq.base(i + 1)][seq.base(j - 1)][seq.base(j - 2)][seq.base(ip)][seq.base(jp)];
            } else if (size1 == 2 and size2 == 1) { /* 1x2 internal loop */
                energy = constants.iloop21[seq.base(jp)][seq.base(ip)][seq.base(j - 1)][seq.base(i + 2)][seq.base(i + 1)][seq.base(j)][seq.base(i)];
            } else if (size == 2) { /* 1*1 internal loops */
                energy = constants.iloop11[seq.base(i)][seq.base(i + 1)][seq.base(ip)][seq.base(j)][seq.base(j - 1)][seq.base(jp)];
            } else if ((size1 == 1 or size2 == 1) and constants.gail) { /* gail = (Grossly Asymmetric Interior Loop Rule) (on/off <-> 1/0)  */
                energy = constants.tstki[seq.base(i)][seq.base(j)][BASE_A][BASE_A] + constants.tstki[seq.base(jp)][seq.base(ip)][BASE_A][BASE_A]
                    + constants.inter[size]
                    + penterm;
            } else { /* General Internal loops */
                energy = constants.tstki[seq.base(i)][seq.base(j)][seq.base(i + 1)][seq.base(j - 1)]
                    + constants.tstki[seq.base(jp)][seq.base(ip)][seq.base(jp + 1)][seq.base(ip - 1)]
                    + constants.inter[size]
                    + penterm;
            }
        }

        return energy;
    }

    Rational NNTM::eH(int i, int j, const RNASequence& seq) const {
        /*
          Compute the energy of a hairpin loop based at pair (i, j)
        */

        // Input specification
        assert (i >= 0 and i < seq.len());
        assert (j >= 0 and j < seq.len());
        assert (i < j);

        Rational energy = Rational::infinity();

        int size = j - i - 1; /*  size is the number of bases in the loop, when the closing pair is excluded */

        /*  look in hairpin, and be careful that there is only 30 values */

        if (size > 30) {
            energy = constants.hairpin[30] + eLL(size) + constants.tstkh[seq.base(i)][seq.base(j)][seq.base(i + 1)][seq.base(j - 1)]; /* size penalty + terminal mismatch stacking energy*/
        }

        else if (size <= 30 and size > 4) {
            energy = constants.hairpin[size] + constants.tstkh[seq.base(i)][seq.base(j)][seq.base(i + 1)][seq.base(j - 1)]; /* size penalty + terminal mismatch stacking energy*/
        }

        else if (size == 4) {
            std::string loopkey = seq.subsequence(i, j);
            Rational tlink = 0; // Loop contribution is typically 0
            if (constants.tloop.count(loopkey) != 0) {
                tlink = constants.tloop.find(loopkey)->second; // But some loops have special contributions, stored in this table
            }
            energy = tlink + constants.hairpin[size] + constants.tstkh[seq.base(i)][seq.base(j)][seq.base(i + 1)][seq.base(j - 1)];
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

        else if (size < 3 and size != 0) {
            /*  no terminal mismatch */
            energy = constants.hairpin[size];
        } else if (size == 0)
            energy = Rational::infinity();

        /*  GGG Bonus => GU closure preceded by GG */
        /*  i-2 = i-1 = i = G, and j = U; i < j */
        if (i >= 2) {
            if (seq.base(i - 2) == BASE_G and seq.base(i - 1) == BASE_G and seq.base(i) == BASE_G
                and seq.base(j) == BASE_U) {
                energy += constants.gubonus;
                /*  printf ("\n GGG bonus for i %d j %d ", i, j); */
            }
        }

        /*  Poly-C loop => How many C are needed for being a poly-C loop */
        Rational tlink = 1;
        for (int index = 1; (index <= size) and (tlink == 1); ++index) {
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

    Rational NNTM::eS(int i, int j, const RNASequence& seq) const {
        /*
          Compute the energy of a stack of pairs (i, j) and (i+1, j-1)
        */
        // Input specification
        assert (i >= 0 and i < seq.len());
        assert (j >= 0 and j < seq.len());
        assert (i < j);

        return constants.stack[seq.base(i)][seq.base(j)][seq.base(i+1)][seq.base(j-1)];
    }

    Rational NNTM::calcVBI(int i, int j, const RNASequenceWithTables& seq) const {
        /*
          Helper method to populate the VBI array
        */

        // Input specification
        assert (i >= 0 and i < seq.len());
        assert (j >= 0 and j < seq.len());
        assert (i < j);

        MinBox<Rational> vals;
        vals.insert(Rational::infinity());

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
                if (q - p > TURN and seq.can_pair(p, q)) {
                    vals.insert(eL(i, j, p, q, seq) + seq.V[p][q]);
                }
            }
        }

        Rational VBIij = vals.minimum();
        return VBIij;
    }
};
