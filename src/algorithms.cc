/*
  GTfold: compute minimum free energy of RNA secondary structure
  Copyright (C) 2008 David A. Bader
  http://www.cc.gatech.edu/~bader
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <assert.h>
#include "constants.h"
#include "utils.h"
#include "energy.h"
#include "global.h"
#include "algorithms.h"
#include "constraints.h"
#include "parametrizer_types.h"
#include <gmpxx.h>
#include <iostream>

#include <vector>
#include <algorithm>

#ifdef _OPENMP
#include "omp.h"
#endif

//#define DEBUG 1

void initializeMatrix(int len) {
    int i, j;

    for (i = 1; i <= len; ++i)
        for (j = len; j >= i; --j)
            if (canPair(RNA[i],RNA[j]) && j-i > TURN)
                PP[i][j]  = 1;
}

void prefilter(int len, int prefilter1, int prefilter2) {
    char** in;
    int i, j, k, count;

    in = (char**)malloc(len*sizeof(char*));
    for (i = 1; i <= len; ++i) in[i - 1] = (char*)malloc(len*sizeof(char));

    for (i = 1; i <= len - prefilter2 + 1; ++i)
        for (j = len; j >= prefilter2 && j >= i; --j) {
            count = 0;
            for (k = 0; k < prefilter2 && k <= (j - i) / 2; ++k)
                if (PP[i + k][j - k] == 1) ++count;
            if (count >= prefilter1)
                for (k = 0; k < prefilter2 && k <= (j - i) / 2; ++k)
                    ++in[i + k - 1][j - k - 1];
        }

    for (i = 1; i <= len; ++i) {
        for (j = len; j >= i; --j)
            if (!in[i - 1][j - 1]) PP[i][j] = 0;
        free(in[i - 1]);
    }

    free(in);
}

mpq_class calcVBI_f(int i, int j) {
    int p=0, q=0;
    mpq_class VBIij = INFINITY_;

    for (p = i+1; p <= MIN(j-2-TURN,i+MAXLOOP+1) ; p++) {
        int minq = j-i+p-MAXLOOP-2;
        if (minq < p+1+TURN) minq = p+1+TURN;
        int maxq = (p==(i+1))?(j-2):(j-1);

        for (q = minq; q <= maxq; q++) {
            if (PP[p][q]==0) continue;
            if (!canILoop(i,j,p,q)) continue;
            std::vector<mpq_class> vals;
            vals.push_back(eL(i, j, p, q) + V_f(p,q));
            vals.push_back(VBIij);
            VBIij = *std::min_element(vals.begin(), vals.end());
        }
    }

    return VBIij;
}

mpq_class calcVBI1(int i, int j) {
    int p=0, q=0;
    mpq_class VBIij = INFINITY_;

    for (p = i+1; p <= MIN(j-2-TURN,i+MAXLOOP+1) ; p++) {
        int minq = j-i+p-MAXLOOP-2;
        if (minq < p+1+TURN) minq = p+1+TURN;
        int maxq = (p==(i+1))?(j-2):(j-1);

        for (q = minq; q <= maxq; q++) {
            if (PP[p][q]==0) continue;
            if (!canILoop(i,j,p,q)) continue;
            std::vector<mpq_class> vals;
            vals.push_back(eL1(i, j, p, q) + V_f(p,q));
            vals.push_back(VBIij);
            VBIij = *std::min_element(vals.begin(), vals.end());
        }
    }

    return VBIij;
}

mpq_class calcVBI2(int i, int j, int  len) {
    int d, ii, jj;
    mpq_class energy = INFINITY_;

    for (d = j-i-3; d >= TURN+1 && d >= j-i-2-MAXLOOP; --d)
        for (ii = i + 1; ii < j - d && ii <= len; ++ii)
        {
            jj = d + ii;
            if (PP[ii][jj]==1) {
                std::vector<mpq_class> vals;
                vals.push_back(energy);
                vals.push_back(eL1(i, j, ii, jj) + V_f(ii, jj));
                energy = *std::min_element(vals.begin(), vals.end());
            }
        }

    return energy;
}

mpq_class calculate(int len) {
    int b, i, j;
#ifdef _OPENMP
    if (g_nthreads > 0) omp_set_num_threads(g_nthreads);
#endif

    initializeMatrix(len);
    if (g_unamode || g_prefilter_mode) {
        prefilter(len,g_prefilter1,g_prefilter2);
    }

    for (b = TURN+1; b <= len-1; b++) {
#ifdef _OPENMP
#pragma omp parallel for private (i,j) schedule(guided)
#endif
        for (i = 1; i <= len - b; i++) {
            j = i + b;

            mpq_class eh, es;
            if (PP[i][j] == 1) {
                if (canHairpin(i,j)) {
                    eh = eH(i,j);
                } else {
                    eh = INFINITY_;
                }

                if (canStack(i,j)) {
                    es = eS(i,j)+V_f(i+1,j-1);
                } else {
                    es = INFINITY_;
                }

                // Internal Loop BEGIN
                if (g_unamode)
                    VBI_f(i,j) = calcVBI1(i,j);
                else
                    VBI_f(i,j) = calcVBI_f(i,j);
                // Internal Loop END

                // Multi Loop BEGIN
                mpq_class d3, d5;
                if (canSS(j-1)) {
                    d3 = Ed3(i,j,j-1);
                } else {
                    d3 = INFINITY_;
                }

                if (canSS(i+1)) {
                    d5 = Ed5(i,j,i+1);
                } else {
                    d5 = INFINITY_;
                }

                if (g_dangles == 2) { // -d2
                    std::vector<mpq_class> vals;
                    vals.push_back(VM_f(i, j));
                    vals.push_back(WMPrime[i+1][j-1] + d3 + d5 + auPenalty(i,j) + multConst[0] + multConst[2]);
                    VM_f(i,j) = *std::min_element(vals.begin(), vals.end());
                } else if (g_dangles == 0) { // -d0
                    std::vector<mpq_class> vals;
                    vals.push_back(VM_f(i,j));
                    vals.push_back(WMPrime[i+1][j-1] + auPenalty(i,j) + multConst[0] + multConst[2]);
                    VM_f(i,j) = *std::min_element(vals.begin(), vals.end());
                } else { // default
                    std::vector<mpq_class> vals;
                    vals.push_back(VM_f(i,j));
                    vals.push_back(WMPrime[i+1][j-1] + auPenalty(i,j) + multConst[0] + multConst[2]);
                    vals.push_back(WMPrime[i+2][j-1] + d5 + auPenalty(i,j) + multConst[0] + multConst[2] + multConst[1]);
                    vals.push_back(WMPrime[i+1][j-2] + d3 + auPenalty(i,j) + multConst[0] + multConst[2] + multConst[1]);
                    vals.push_back(WMPrime[i+2][j-2] + d3 + d5 + auPenalty(i,j) + multConst[0] + multConst[2] + 2*multConst[1]);
                    VM_f(i,j) = *std::min_element(vals.begin(), vals.end());
                }

                if (!canStack(i,j)) VM_f(i,j) = INFINITY_;

                // Multi Loop END

                std::vector<mpq_class> vals;
                vals.push_back(eh);
                vals.push_back(es);
                vals.push_back(VBI_f(i,j));
                vals.push_back(VM_f(i,j));
                V_f(i,j) = *std::min_element(vals.begin(), vals.end());
            }
            else {
                V_f(i,j) = INFINITY_;
            }

            // Added auxillary storage WMPrime to speedup multiloop calculations
            int h;
            for (h = i+TURN+1 ; h <= j-TURN-2; h++) {
                std::vector<mpq_class> vals;
                vals.push_back(WMPrime[i][j]);
                vals.push_back(WMU_f(i,h-1) + WML_f(h,j));
                WMPrime[i][j] = *std::min_element(vals.begin(), vals.end());

            }

            // WM begin
            mpq_class newWM = INFINITY_;

            //ZS: This sum corresponds to when i,j are NOT paired with each other.
            //So we need to make sure only terms where i,j aren't pairing are considered.
            if (!forcePair(i,j)) {
                std::vector<mpq_class> vals;
                vals.push_back(newWM);
                vals.push_back(WMPrime[i][j]);
                newWM = *std::min_element(vals.begin(), vals.end());
            }

            if (g_dangles == 2) {
                mpq_class energy = V_f(i,j) + auPenalty(i,j) + multConst[2];
                if (i == 1) {
                    energy += Ed3(j,i,len);
                } else {
                    energy += Ed3(j,i,i-1);
                }

                /*if (j<len)*/ energy += Ed5(j,i,j+1);
                if (canSS(i)&&canSS(j)){
                    std::vector<mpq_class> vals;
                    vals.push_back(energy);
                    vals.push_back(newWM);
                    newWM = *std::min_element(vals.begin(), vals.end());
                }
            } else if (g_dangles == 0) {
                std::vector<mpq_class> vals;
                vals.push_back(V_f(i,j) + auPenalty(i,j) + multConst[2]);
                vals.push_back(newWM);
                newWM = *std::min_element(vals.begin(), vals.end());
            } else { // default
                std::vector<mpq_class> vals;
                vals.push_back(newWM);
                vals.push_back(V_f(i,j) + auPenalty(i,j) + multConst[2]);

                if (canSS(i))
                    vals.push_back(V_f(i+1,j) + Ed3(j,i+1,i) + auPenalty(i+1,j) + multConst[2] + multConst[1]); //i dangle

                if (canSS(j))
                    vals.push_back(V_f(i,j-1) + Ed5(j-1,i,j) + auPenalty(i,j-1) + multConst[2] + multConst[1]);  //j dangle

                if (canSS(i)&&canSS(j))
                    vals.push_back(V_f(i+1,j-1) + Ed3(j-1,i+1,i) + Ed5(j-1,i+1,j) + auPenalty(i+1,j-1) + multConst[2] + 2*multConst[1]); //i,j dangle

                newWM = *std::min_element(vals.begin(), vals.end());
            }

            std::vector<mpq_class> vals;
            vals.push_back(newWM);

            if (canSS(i))
                vals.push_back(WMU_f(i+1,j) + multConst[1]); //i dangle

            if (canSS(j))
                vals.push_back(WML_f(i,j-1) + multConst[1]); //j dangle

            newWM = *std::min_element(vals.begin(), vals.end());

            WMU_f(i,j) = WML_f(i,j) = newWM;
            // WM end
        }
    }

    W[0] = 0;
    for (j = 1; j <= len; j++) {
        int i;
        mpq_class Wj, Widjd, Wijd, Widj, Wij, Wim1;
        Wj = 0;
        for (i = 1; i < j-TURN; i++) {
            Wij = Widjd = Wijd = Widj = INFINITY_;
            Wim1 = MIN(0, W[i-1]);

            if (g_dangles == 2) { // -d2 option
                mpq_class energy = V_f(i,j) +	 auPenalty(i,j) + Wim1;
                if (i>1) energy +=  Ed3(j,i,i-1);
                if (j<len) energy += Ed5(j,i,j+1);
                if (canSS(i)&&canSS(j)) Widjd = energy;

                std::vector<mpq_class> vals;
                vals.push_back(Wij);
                vals.push_back(Widjd);

                Wij = *std::min_element(vals.begin(), vals.end());
            }	else if (g_dangles == 0) { // -d0 option
                Wij = V_f(i, j) + auPenalty(i, j) + Wim1;
            } else { // default
                Wij = V_f(i, j) + auPenalty(i, j) + Wim1;
                if (canSS(i)) {
                    Widj = V_f(i+1, j) + auPenalty(i+1,j) + Ed3(j,i + 1,i) + Wim1;
                } else {
                    Widj = INFINITY_;
                }

                if (canSS(j)) {
                    Wijd = V_f(i,j-1) + auPenalty(i,j-1) + Ed5(j-1,i,j) + Wim1;
                } else {
                    Wijd = INFINITY_;
                }

                if (canSS(i)&&canSS(j)) {
                    Widjd = V_f(i+1,j-1) + auPenalty(i+1,j-1) + Ed3(j-1,i + 1,i) + Ed5(j-1,i+1,j) + Wim1;
                } else {
                    Widjd = INFINITY_;
                }

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
        if (canSS(j))
            vals.push_back(W[j-1]);

        W[j] = *std::min_element(vals.begin(), vals.end());
    }

#ifdef DEBUG
    FILE* file = fopen("VBI.txt", "w");
    int ii, jj;
    for (ii = 1; ii <= len; ++ii) {
        for (jj = len; jj > ii; --jj) {
            fprintf(file, "%d %d %f\n",ii,jj,VBI_f(ii,jj).get_d());
        }
    }
    fclose(file);

    file = fopen("Eh.txt", "w");
    for (ii = 1; ii <= len; ++ii) {
        for (jj = len; jj > ii; --jj) {
            mpq_class eh = INFINITY_;
            if (PP[ii][jj])	eh = eH(ii,jj);
            fprintf(file, "%d %d %f\n",ii,jj,eh>=INFINITY_?INFINITY_.get_d():eh.get_d());
        }
    }
    fclose(file);

    file = fopen("Es.txt", "w");
    for (ii = 1; ii <= len; ++ii) {
        for (jj = len; jj > ii; --jj) {
            mpq_class es = INFINITY_;
            if (PP[ii][jj] && PP[ii+1][jj-1]) es = eS(ii,jj);
            fprintf(file, "%d %d %f\n",ii,jj,es>=INFINITY_?INFINITY_.get_d():es.get_d());
        }
    }
    fclose(file);

    file = fopen("BP.txt", "w");
    for (ii = 1; ii <= len; ++ii) {
        for (jj = len; jj > ii; --jj) {
            fprintf(file, "%d %d %d\n",ii,jj,PP[ii][jj]);
        }
    }
    fclose(file);

    file = fopen("VM.txt", "w");
    for (ii = 1; ii <= len; ++ii) {
        for (jj = len; jj > ii; --jj) {
            fprintf(file, "%d %d %f\n",ii,jj,VM_f(ii,jj).get_d());
        }
    }
    fclose(file);

    file = fopen("WM.txt", "w");
    for (ii = 1; ii <= len; ++ii) {
        for (jj = len; jj > ii; --jj) {
            fprintf(file, "%d %d %f\n",ii,jj,WM_f(ii,jj).get_d());
        }
    }
    fclose(file);
#endif

    return W[len];
}
