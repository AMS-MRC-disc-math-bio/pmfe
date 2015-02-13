/**
   GTfold: compute minimum free energy of RNA secondary structure
   Copyright (C) 2008  David A. Bader
   http://www.cc.gatech.edu/~bader

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   @author prashant {pgaurav@gatech.edu}

*/

#include <cstdio>

#include "constants.h"
#include "energy.h"
#include "utils.h"
#include "global.h"
#include "subopt_traceback.h"

#include <iostream>
#include <iterator>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <gmpxx.h>

namespace pmfe {
    using std::pair;
    using std::cout;
    using std::endl;

    const char* lstr[] = {"W", "V", "VBI", "VM", "WM", "WMPrime", "M", "M1"};
    const char d3symb = '>';
    const char d5symb = '<';

    void (*trace_func[8]) (int i, int j, ps_t& ps, ps_stack_t& gs);
    static mpq_class mfe_ = INFINITY_;
    static mpq_class delta_ = 0;
    static int length = -1;
    static int gflag = 0;

    bool DEBUG = false;

    //TODO: Make these dynamic to eliminate hard ceiling at len=1500
    static mpq_class FM1[1500][1500] = {{mpq_class(0)}};
    static mpq_class FM[1500][1500] = {{mpq_class(0)}};

    // k dangles on 5' end of (i, j)
    static inline mpq_class Ed5_new(int i, int j, int k) {
        return (k!=0) ? Ed3(j, i, k) : Ed3(j,i,length);
    }

    // k dangles on 3' end of (i, j)
    static inline mpq_class Ed3_new(int i, int j, int k) {
        return (k!=length+1)?Ed5(j, i, k) : Ed5(j, i, 1);
    }

    // dangle on the 5' end of (i, j)
    // if inside==true, dangle i+1 instead of i-1
    static inline mpq_class Ed5_pair(int i, int j, bool inside = false) {
        if (i >= 2) {
            if (inside)
                return Ed5(i, j, i+1);
            else
                return Ed3(j, i, i-1);
        } else if (i == 1) {
            if (inside)
                return Ed5(i, j, i+1);
            else
                return Ed3(j, i, length);
        } else {
            exit(EXIT_FAILURE);
        }
    }

    // dangle on the 3' end of (i, j)
    // if inside==true, dangle j-1 instead of j+1
    static inline mpq_class Ed3_pair(int i, int j, bool inside = false) {
        if (j <= length-1) {
            if (inside)
                return Ed3(i, j, j-1);
            else
                return Ed5(j, i, j+1);
        } else if (j == length) {
            if (inside)
                return Ed3(i, j, j-1);
            else
                return Ed5(j, i, 1);
        } else {
            exit(EXIT_FAILURE);
        }
    }

    // Per Wuchty et al., compute MFE for a final branch in a multiloop which
    // begins at i and ends at j, possibly including free bases at the 5' end
    void calculate_fm1() {

        for (int i = 1; i <= length; ++i) {
            for (int j = i+1; j <= length; ++j) {
                mpq_class min = INFINITY_;
                for (int l = i+TURN+1; l <= j; ++l) {

                    mpq_class fm1 = 0;
                    switch (g_dangles) {
                    case NO_DANGLE:
                        {
                            fm1 = V_f(i,l) + auPenalty(i,l) + multConst[1]*(j-l) + multConst[2];
                            break;
                        }

                    case CHOOSE_DANGLE:
                        {
                            // Experimental d1 version; allows the stem start to move so the dangles are included in the substructure
                            mpq_class d5 = Ed5_pair(i+1,l);
                            mpq_class d3 = Ed3_pair(i,l-1);
                            mpq_class d53 = Ed5_pair(i+1, l-1) + Ed3_pair(i+1, l-1);
                            std::vector<mpq_class> vals;
                            vals.push_back(V_f(i,l) + auPenalty(i,l) + multConst[1]*(j - l) + multConst[2]);
                            vals.push_back(V_f(i+1,l) + auPenalty(i+1,l) + d5 + multConst[1]*(j - l + 1) + multConst[2]);
                            vals.push_back(V_f(i,l-1) + auPenalty(i,l-1) + d3 + multConst[1]*(j - (l-1)) + multConst[2]);
                            vals.push_back(V_f(i+1,l-1) + auPenalty(i+1,l-1) + d53 + multConst[1]*(j - (l-1) + 1) + multConst[2]);
                            fm1 = *std::min_element(vals.begin(), vals.end());
                            break;
                        }

                    case BOTH_DANGLE:
                        {
                            mpq_class d5 = Ed5_new(i,l,i-1);
                            mpq_class d3 = Ed3_new(i,l,l+1);
                            fm1 = V_f(i,l) + auPenalty(i,l) + d5 + d3 + multConst[1]*(j-l) + multConst[2];
                            break;
                        }

                    default:
                        {
                            exit(EXIT_FAILURE);
                            break;
                        }
                    }
                    min = std::min(min, fm1);
                }
                FM1[i][j] = min;
            }
        }

    }

    void calculate_fm() {
        for (int i = 1; i <= length; ++i) {
            for (int j = i+1; j <= length;++j) {
                mpq_class min1 = INFINITY_;
                for (int k = i+TURN+1; k <= j-TURN-1; ++k) {
                    mpq_class x = FM[i][k-1] + FM1[k][j];
                    min1 = std::min(min1, x);
                }
                mpq_class min2 = INFINITY_;
                for (int k = i; k <= j-TURN-1; ++k) {
                    mpq_class x = FM1[k][j] + multConst[1]*(k-i);
                    min2 = std::min(min2, x);
                }
                FM[i][j] = std::min(min1, min2);
            }
        }
    }

    void process(ss_map_t& subopt_data, int len, string suboptFile, int max_structure_count) {
        ofstream outfile;
        outfile.open(suboptFile.c_str(), ios::out | ios::app);
        char buff[4096]; // TODO: make dynamic to avoid hard ceiling

        int count = 0;
        length = len;

        calculate_fm1();
        calculate_fm();

        ps_stack_t gstack;

        // initialize the partial structure, segment stack = {[1,n]}, label = W, list_bp = {}
        ps_t first(0, len);
        first.push(segment(1, len, lW, W[len]));
        gstack.push(first); // initialize the partial structure stack

        // TODO: parallelize?
        while (!gstack.empty()) {
            ps_t ps = gstack.top();
            gstack.pop();

            if (DEBUG) printf("Starting loop with sequence %s\n", ps.str.c_str());

            if (ps.empty()) {
                count++;
                if (DEBUG) printf("Sequence complete; writing!\nOutput:%d\t%s\t%f\n", count, (ps.str).c_str(), ps.ae_.get_d());
                //cout << ps.str << endl;
                sprintf(buff,"%d\t%s\t%f", count, (ps.str).c_str(), ps.ae_.get_d());
                outfile << buff << std::endl;
                if(max_structure_count>0 && count>=max_structure_count) break;//exit
                continue;
            }
            else {
                segment smt = ps.top();
                ps.pop();

                gflag = 0;
                if (smt.j_ - smt.i_ > TURN) {
                    if (DEBUG) printf("Segment (%i, %i), trace %s\n", smt.i_, smt.j_, lstr[smt.label_]);
                    (*trace_func[smt.label_])(smt.i_, smt.j_, ps, gstack);
                } else {
                    if (DEBUG) printf("Segment (%i, %i) too short for trace %s\n", smt.i_, smt.j_, lstr[smt.label_]);
                }

                // discarded current segment, using remaining ones
                if (!gflag) {
                    if (DEBUG) printf("Completed segment (%i, %i)\n", smt.i_, smt.j_);
                    ps_t ps1(ps);
                    gstack.push(ps1);
                }
            }
        }
        outfile.close();
        printf("Counts of structure generated=%d\n", count);

#ifdef DEBUG
        if (DEBUG) printf("# SS = %d\n", count);
#endif
    }

    ss_map_t subopt_traceback(int len, mpq_class delta, string suboptFile,  int max_structure_count) {
        trace_func[lW] = traceW;
        trace_func[lV] = traceV;
        trace_func[lVBI] = traceVBI;
        //trace_func[lVM] = traceVM;
        //trace_func[lWM] = traceWM;
        //trace_func[lWMPrime] = traceWMPrime;
        trace_func[lM] = traceM;
        trace_func[lM1] = traceM1;

        mfe_ = W[len];
        delta_ = delta;

        ss_map_t subopt_data;
        process(subopt_data, len, suboptFile, max_structure_count);

        return subopt_data;
    }

    void traceV(int i, int j, ps_t& ps, ps_stack_t& gstack) {
        if (DEBUG) printf("traceV(%i, %i):\n\tpartial structure %s\n", i, j, ps.str.c_str());
        // Hairpin Loop
        if (eH(i,j) + ps.total()  <= mfe_ + delta_) {
            if (DEBUG) printf("Hairpin (%i, %i)\n", i, j);
            ps_t ps1(ps);
            ps1.accumulate(eH(i,j));
            ps1.update(i, j, '(', ')');
            push_to_gstack(gstack, ps1);
        }

        // Stack
        if (eS(i, j) + V_f(i+1, j-1) + ps.total() <= mfe_ + delta_) {
            if (DEBUG) printf("Stack (%i, %i)\n", i, j);
            ps_t ps1(ps);
            ps1.push(segment(i+1, j-1, lV, V_f(i+1, j-1)));
            ps1.accumulate(eS(i,j));
            ps1.update(i, j , '(', ')');
            push_to_gstack(gstack, ps1);
        }

        // Internal Loop
        if (VBI_f(i,j) + ps.total() <= mfe_ + delta_) {
            if (DEBUG) printf("Intloop (%i, %i)\n", i, j);
            traceVBI(i,j,ps,gstack);
        }

        // Multiloop

        int k;

        for (k = i+2; k <= j-TURN-1; ++k) {
            switch (g_dangles) {
            case NO_DANGLE:
                {
                    mpq_class kenergy1 = FM[i+1][k] + FM1[k+1][j-1];
                    mpq_class kenergy2 = auPenalty(i,j) + multConst[0] + multConst[2];
                    mpq_class kenergy_total = kenergy1 + kenergy2;
                    if (DEBUG) printf("Multiloop check (%i, %i, %i): %8.2f <= %8.2f\n", i, k, j, kenergy_total.get_d() + ps.total().get_d(), mfe_.get_d() + delta_.get_d());
                    if (kenergy_total + ps.total() <= mfe_ + delta_) {
                        if (DEBUG) printf("Multiloop confirmed at (%i, %i, %i)\n", i, k, j);
                        ps_t ps1(ps);
                        ps1.push(segment(i+1,k, lM, FM[i+1][k]));
                        ps1.push(segment(k+1,j-1, lM1, FM1[k+1][j-1]));
                        ps1.accumulate(kenergy2);
                        ps1.update(i,j,'(',')');
                        push_to_gstack(gstack, ps1);
                    }
                    break;
                }

            case CHOOSE_DANGLE:
                {
                    mpq_class d5 = Ed5_pair(i, j, true);
                    mpq_class d3 = Ed3_pair(i, j, true);
                    mpq_class d53 = d5 + d3;
                    if (FM[i+1][k] + FM1[k+1][j-1] + auPenalty(i,j) + multConst[0] + multConst[2] + ps.total() <= mfe_ + delta_) {
                        printf("traceV multiloop no dangle (%i, %i, %i)\n", i, k, j);
                        ps_t ps1(ps);
                        ps1.push(segment(i+1, k, lM, FM[i+1][k]));
                        ps1.push(segment(k+1, j-1, lM1, FM1[k+1][j-1]));
                        ps1.accumulate(auPenalty(i,j) + multConst[0] + multConst[2]);
                        ps1.update(i,j,'(',')');
                        push_to_gstack(gstack, ps1);
                    }
                    if (FM[i+2][k] + FM1[k+1][j-1] + auPenalty(i,j) + d5 + multConst[0] + multConst[1] + multConst[2] + ps.total() <= mfe_ + delta_) {
                        printf("traceV multiloop d5 (%i, %i, %i)\n", i, k, j);
                        ps_t ps1(ps);
                        ps1.push(segment(i+2,k, lM, FM[i+2][k]));
                        ps1.push(segment(k+1,j-1, lM1, FM1[k+1][j-1]));
                        ps1.accumulate(auPenalty(i,j) + d5 + multConst[0] + multConst[1] + multConst[2]);
                        ps1.update(i,j,'(',')');
                        ps1.update(i+1,d3symb);
                        push_to_gstack(gstack, ps1);
                    }
                    if (FM[i+1][k] + FM1[k+1][j-2] + auPenalty(i,j) + d3 + multConst[0] + multConst[1] + multConst[2] + ps.total() <= mfe_ + delta_) {
                        printf("traceV multiloop d3 (%i, %i, %i)\n", i, k, j);
                        ps_t ps1(ps);
                        ps1.push(segment(i+1,k, lM, FM[i+1][k]));
                        ps1.push(segment(k+1,j-2, lM1, FM1[k+1][j-2]));
                        ps1.accumulate(auPenalty(i,j) + d3 + multConst[0] + multConst[1] + multConst[2]);
                        ps1.update(i,j,'(',')');
                        ps1.update(j,d5symb);
                        push_to_gstack(gstack, ps1);
                    }
                    if (FM[i+2][k] + FM1[k+1][j-2] + auPenalty(i,j) + d53 + multConst[0] + 2*multConst[1] + multConst[2] + ps.total() <= mfe_ + delta_) {
                        printf("traceV multiloop d53 (%i, %i, %i)\n", i, k, j);
                        ps_t ps1(ps);
                        ps1.push(segment(i+2,k, lM, FM[i+2][k]));
                        ps1.push(segment(k+1,j-2, lM1, FM1[k+1][j-2]));
                        ps1.accumulate(auPenalty(i,j) + d53 + multConst[0] + 2*multConst[1] + multConst[2]);
                        ps1.update(i,j,'(',')');
                        ps1.update(i,j,d3symb,d5symb);
                        push_to_gstack(gstack, ps1);
                    }
                    break;
                }

            case BOTH_DANGLE:
                {
                    mpq_class d5 = Ed5(i, j, i+1);
                    mpq_class d3 = Ed3(i, j, j-1);
                    mpq_class kenergy1 = FM[i+1][k] + FM1[k+1][j-1];
                    mpq_class kenergy2 = d5 + d3 + auPenalty(i,j) + multConst[0] + multConst[2];
                    mpq_class kenergy_total = kenergy1 + kenergy2;
                    if (kenergy_total + ps.total() <= mfe_ + delta_) {
                        if (DEBUG) printf("traceV multiloop (%i, %i, %i)\n", i, k, j);
                        ps_t ps1(ps);
                        ps1.push(segment(i+1,k, lM, FM[i+1][k]));
                        ps1.push(segment(k+1,j-1, lM1, FM1[k+1][j-1]));
                        ps1.accumulate(kenergy2);
                        ps1.update(i,j,'(',')');
                        push_to_gstack(gstack, ps1);
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
    }

    void traceVBI(int i, int j, ps_t& ps, ps_stack_t& gstack) {
        if (DEBUG) printf("traceVBI(%i, %i):\n\tpartial structure %s\n", i, j, ps.str.c_str());
        int p,q;

        for (p = i+1; p <= std::min(j-2-TURN,i+MAXLOOP+1) ; p++) {
            int minq = j-i+p-MAXLOOP-2;
            if (minq < p+1+TURN) minq = p+1+TURN;
            int maxq = (p==(i+1))?(j-2):(j-1);
            for (q = minq; q <= maxq; q++) {
                if (V_f(p, q) + eL(i, j, p, q) + ps.total() <= mfe_ + delta_) {
                    ps_t ps1(ps);
                    ps1.push(segment(p, q, lV, V_f(p, q)));
                    ps1.update(i, j , '(', ')');
                    ps1.accumulate(eL(i, j, p, q));
                    push_to_gstack(gstack, ps1);
                }
            }
        }
    }

    void traceW(int i, int j, ps_t& ps, ps_stack_t& gstack) {
        if (DEBUG) printf("traceW(%i, %i):\n\tpartial structure %s\n", i, j, ps.str.c_str());
        for (int l = i; l < j-TURN; ++l) {
            mpq_class wim1 = W[l-1];

            switch (g_dangles){
            case NO_DANGLE:
                {
                    mpq_class bonus = auPenalty(l, j);
                    if (V_f(l,j) + wim1 + bonus + ps.total() <= mfe_ + delta_ ) {
                        if (DEBUG) printf("traceW: (%i, %i, %i)\n", i, l, j);
                        ps_t ps1(ps);
                        ps1.push(segment(l, j, lV, V_f(l,j)));
                        if (l > i) ps1.push(segment(i, l-1, lW, wim1));
                        ps1.accumulate(bonus);
                        push_to_gstack(gstack, ps1);
                    }
                    break;
                }

            case CHOOSE_DANGLE:
                {
                    mpq_class d5 = Ed5_pair(l+1, j);
                    mpq_class d3 = Ed3_pair(l, j-1);
                    mpq_class d53 = Ed5_pair(l+1, j-1) + Ed3_pair(l+1, j-1);

                    if (V_f(l,j) + auPenalty(l, j) + wim1 + ps.total() <= mfe_ + delta_) {
                        if (DEBUG) printf("traceW no dangle: (%i, %i, %i)\n", i, l, j);
                        ps_t ps1(ps);
                        ps1.push(segment(l, j, lV, V_f(l,j)));
                        if (l > i) ps1.push(segment(i, l-1, lW, wim1));
                        ps1.accumulate(auPenalty(l, j));
                        push_to_gstack(gstack, ps1);
                    }
                    if (V_f(l+1,j) + auPenalty(l+1,j) + multConst[1] + d5 + wim1 + ps.total() <= mfe_ + delta_) {
                        if (DEBUG) printf("traceW d5: (%i, %i, %i)\n", i, l, j);
                        ps_t ps1(ps);
                        ps1.push(segment(l+1, j, lV, V_f(l+1,j)));
                        ps1.update(l, d5symb);
                        if (l > i) ps1.push(segment(i, l-1, lW, wim1));
                        ps1.accumulate(auPenalty(l+1, j) + multConst[1] + d5);
                        push_to_gstack(gstack, ps1);
                    }
                    if (V_f(l,j-1) + auPenalty(l,j-1) + multConst[1] + d3 + wim1 + ps.total() <= mfe_ + delta_) {
                        if (DEBUG) printf("traceW d3: (%i, %i, %i)\n", i, l, j);
                        ps_t ps1(ps);
                        ps1.push(segment(l, j-1, lV, V_f(l,j-1)));
                        ps1.update(j, d3symb);
                        if (l > i) ps1.push(segment(i, l-1, lW, wim1));
                        ps1.accumulate(auPenalty(l, j-1) + multConst[1] + d3);
                        push_to_gstack(gstack, ps1);
                    }
                    if (V_f(l+1,j-1) + auPenalty(l+1, j-1) + 2*multConst[1] + d53 + wim1 + ps.total() <= mfe_ + delta_) {
                        if (DEBUG) printf("traceW d53: (%i, %i, %i)\n", i, l, j);
                        ps_t ps1(ps);
                        ps1.push(segment(l+1, j-1, lV, V_f(l+1,j-1)));
                        ps1.update(l, j, d5symb, d3symb);
                        if (l > i) ps1.push(segment(i, l-1, lW, wim1));
                        ps1.accumulate(auPenalty(l+1, j-1) + 2*multConst[1] + d53);
                        push_to_gstack(gstack, ps1);
                    }
                    break;
                }

            case BOTH_DANGLE:
                {
                    mpq_class d3 = (l>i)?Ed3(j,l,l-1):mpq_class(0);
                    mpq_class d5 = (j<length)?Ed5(j,l,j+1):mpq_class(0);
                    mpq_class bonus = auPenalty(l, j) + d5 + d3;
                    if (V_f(l,j) + wim1 + bonus + ps.total() <= mfe_ + delta_ ) {
                        ps_t ps1(ps);
                        ps1.push(segment(l, j, lV, V_f(l,j)));
                        if (l > i) ps1.push(segment(i, l-1, lW, wim1));
                        ps1.accumulate(bonus);
                        push_to_gstack(gstack, ps1);
                    }
                    break;
                }

            default:
                exit(EXIT_FAILURE);
                break;
            };
        };

        if (W[j-1] + ps.total() <= mfe_ + delta_) {
            if (DEBUG) printf("traceW final branch: (%i, %i)\n", i, j);
            ps_t ps1(ps);
            ps1.push(segment(i, j-1, lW, W[j-1]));
            push_to_gstack(gstack, ps1);
        }
    }

    void traceM1(int i, int j, ps_t& ps, ps_stack_t& gstack) {
        if (DEBUG) printf("traceM1(%i, %i):\n\tpartial structure %s\n", i, j, ps.str.c_str());
        if (FM1[i][j-1] + multConst[1] + ps.total() <= mfe_ + delta_) {
            if (DEBUG) printf("traceM1 nibble: (%i, %i)\n", i, j);
            ps_t ps1(ps);
            ps1.push(segment(i, j-1, lM1, FM1[i][j-1]));
            ps1.accumulate(multConst[1]);
            push_to_gstack(gstack, ps1);
        }

        mpq_class bonus = 0;

        switch (g_dangles) {
        case NO_DANGLE:
            {
                mpq_class aup = auPenalty(i,j);
                bonus = aup + multConst[2];
                if (V_f(i,j) + bonus + ps.total() <= mfe_ + delta_) {
                    ps_t ps1(ps);
                    ps1.push(segment(i, j, lV, V_f(i,j)));
                    ps1.accumulate(bonus);
                    ps1.update(i,j,'(',')');
                    push_to_gstack(gstack, ps1);
                }
                break;
            }

        case CHOOSE_DANGLE:
            {
                mpq_class d5 = Ed5_pair(i+1, j);
                mpq_class d3 = Ed3_pair(i, j-1);
                mpq_class d53 = Ed5_pair(i+1, j-1) + Ed3_pair(i+1, j-1);
                if (V_f(i,j) + auPenalty(i,j) + multConst[2] + ps.total() <= mfe_ + delta_) {
                    if (DEBUG) printf("traceM1 no dangle: (%i, %i)\n", i, j);
                    ps_t ps1(ps);
                    ps1.push(segment(i, j, lV, V_f(i,j)));
                    ps1.accumulate(auPenalty(i,j) + multConst[2]);
                    ps1.update(i,j,'(',')');
                    push_to_gstack(gstack, ps1);
                }
                if (V_f(i+1,j) + auPenalty(i+1,j) + multConst[2] + multConst[1] + d5 + ps.total() <= mfe_ + delta_) {
                    if (DEBUG) printf("traceM1 d5: (%i, %i)\n", i, j);
                    ps_t ps1(ps);
                    ps1.push(segment(i+1, j, lV, V_f(i+1,j)));
                    ps1.accumulate(auPenalty(i+1,j) + multConst[2] + multConst[1] + d5);
                    ps1.update(i+1,j,'(',')');
                    ps1.update(i, d5symb);
                    push_to_gstack(gstack, ps1);
                }
                if (V_f(i,j-1) + auPenalty(i,j-1) + multConst[2] + multConst[1] + d3 + ps.total() <= mfe_ + delta_) {
                    if (DEBUG) printf("traceM1 d3: (%i, %i)\n", i, j);
                    ps_t ps1(ps);
                    ps1.push(segment(i, j-1, lV, V_f(i,j-1)));
                    ps1.accumulate(auPenalty(i,j-1) + multConst[2] + multConst[1] + d3);
                    ps1.update(i,j-1,'(',')');
                    ps1.update(j, d3symb);
                    push_to_gstack(gstack, ps1);
                }
                if (V_f(i+1,j-1) + auPenalty(i+1,j-1) + multConst[2] + 2*multConst[1] + d53 + ps.total() <= mfe_ + delta_) {
                    if (DEBUG) printf("traceM1 d53: (%i, %i)\n", i, j);
                    ps_t ps1(ps);
                    ps1.push(segment(i+1, j-1, lV, V_f(i+1,j-1)));
                    ps1.accumulate(auPenalty(i+1,j-1) + multConst[2] + 2*multConst[1] + d53);
                    ps1.update(i+1,j-1,'(',')');
                    ps1.update(i, j, d5symb, d3symb);
                    push_to_gstack(gstack, ps1);
                }
                break;
            }

            //case CHOOSE_DANGLE:
        case BOTH_DANGLE:
            {
                mpq_class d5 = Ed5_new(i, j, i-1);
                mpq_class d3 = Ed3_new(i, j, j+1);
                mpq_class aup = auPenalty(i,j);
                bonus = d5 + d3 + aup + multConst[2];
                if (V_f(i,j) + bonus + ps.total() <= mfe_ + delta_) {
                    ps_t ps1(ps);
                    ps1.push(segment(i, j, lV, V_f(i,j)));
                    ps1.accumulate(bonus);
                    ps1.update(i,j,'(',')');
                    push_to_gstack(gstack, ps1);
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

    void traceM(int i, int j, ps_t& ps, ps_stack_t& gstack) {
        if (DEBUG) printf("traceM(%i, %i):\n\tpartial structure %s\n", i, j, ps.str.c_str());

        if (FM[i][j-1] + multConst[1] + ps.total() <= mfe_ + delta_) {
            if (DEBUG) printf("Multiloop (%i, %i)\n", i, j);
            ps_t ps1(ps);
            ps1.push(segment(i, j-1, lM, FM[i][j-1]));
            ps1.accumulate(multConst[1]);
            push_to_gstack(gstack, ps1);
        }

        mpq_class bonus = 0;

        // case that this whole region is a single branch
        switch (g_dangles) {
        case NO_DANGLE:
            {
                mpq_class aup = auPenalty (i,j);
                bonus = multConst [2] + aup;
                if (V_f(i,j) + bonus + ps.total() <= mfe_ + delta_) {
                    ps_t ps1(ps);
                    ps1.push(segment(i, j, lV, V_f(i,j)));
                    ps1.accumulate(bonus);
                    ps1.update(i,j,'(',')');
                    push_to_gstack(gstack, ps1);
                }
                break;
            }

        case CHOOSE_DANGLE:
            {
                mpq_class d5 = Ed5_pair(i+1, j);
                mpq_class d3 = Ed3_pair(i, j-1);
                mpq_class d53 = Ed5_pair(i+1, j-1) + Ed3_pair(i+1, j-1);
                if (V_f(i,j) + multConst[2] + auPenalty(i,j) + ps.total() <= mfe_ + delta_) {
                    ps_t ps1(ps);
                    ps1.push(segment(i, j, lV, V_f(i,j)));
                    ps1.accumulate(multConst[2] + auPenalty(i,j));
                    ps1.update(i,j,'(',')');
                    push_to_gstack(gstack, ps1);
                }
                if (V_f(i+1,j) + multConst[2] + multConst[1] + auPenalty(i+1,j) + d5 + ps.total() <= mfe_ + delta_) {
                    ps_t ps1(ps);
                    ps1.push(segment(i+1, j, lV, V_f(i+1,j)));
                    ps1.accumulate(multConst[2] + multConst[1] + auPenalty(i+1,j) + d5);
                    ps1.update(i+1,j,'(',')');
                    ps1.update(i, d5symb);
                    push_to_gstack(gstack, ps1);
                }
                if (V_f(i,j-1) + multConst[2] + multConst[1] + auPenalty(i,j-1) + d3 + ps.total() <= mfe_ + delta_) {
                    ps_t ps1(ps);
                    ps1.push(segment(i, j-1, lV, V_f(i,j-1)));
                    ps1.accumulate(multConst[2] + multConst[1] + auPenalty(i,j-1) + d3);
                    ps1.update(i,j-1,'(',')');
                    ps1.update(j, d3symb);
                    push_to_gstack(gstack, ps1);
                }
                if (V_f(i+1,j-1) + multConst[2] + 2*multConst[1] + auPenalty(i+1,j-1) + d53 + ps.total() <= mfe_ + delta_) {
                    ps_t ps1(ps);
                    ps1.push(segment(i+1, j-1, lV, V_f(i+1,j-1)));
                    ps1.accumulate(multConst[2] + 2*multConst[1] + auPenalty(i+1,j-1) + d53);
                    ps1.update(i+1,j-1,'(',')');
                    ps1.update(i, j, d5symb, d3symb);
                    push_to_gstack(gstack, ps1);
                }
                break;
            }

            //case CHOOSE_DANGLE:
        case BOTH_DANGLE:
            {
                mpq_class d5 = Ed5_new (i, j, i-1);
                mpq_class d3 = Ed3_new (i, j, j+1);
                mpq_class aup = auPenalty (i,j);
                bonus = d5 + d3 + multConst [2] + aup;
                if (V_f(i,j) + bonus + ps.total() <= mfe_ + delta_) {
                    ps_t ps1(ps);
                    ps1.push(segment(i, j, lV, V_f(i,j)));
                    ps1.accumulate(bonus);
                    ps1.update(i,j,'(',')');
                    push_to_gstack(gstack, ps1);
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

            mpq_class bonus = 0;

            switch (g_dangles) {
            case NO_DANGLE:
                {
                    mpq_class aup = auPenalty(k+1, j);
                    bonus = multConst[2] + aup;
                    if (FM[i][k] + V_f(k+1,j) + bonus + ps.total() <= mfe_ + delta_) {
                        ps_t ps1(ps);
                        ps1.push(segment(i, k, lM, FM[i][k]));
                        ps1.push(segment(k+1, j, lV, V_f(k+1,j)));
                        ps1.accumulate(bonus);
                        ps1.update(k+1,j,'(',')');
                        push_to_gstack(gstack, ps1);
                    }
                    break;
                }

            case CHOOSE_DANGLE:
                {
                    mpq_class d5 = Ed5_pair(k+2, j);
                    mpq_class d3 = Ed3_pair(k+1, j-1);
                    mpq_class d53 = Ed5_pair(k+2, j-1) + Ed3_pair(k+2, j-1);
                    if (FM[i][k] + V_f(k+1,j) + multConst[2] + auPenalty(k+1, j) + ps.total() <= mfe_ + delta_) {
                        ps_t ps1(ps);
                        ps1.push(segment(i, k, lM, FM[i][k]));
                        ps1.push(segment(k+1, j, lV, V_f(k+1,j)));
                        ps1.accumulate(multConst[2] + auPenalty(k+1, j));
                        ps1.update(k+1,j,'(',')');
                        push_to_gstack(gstack, ps1);
                    }
                    if (FM[i][k] + V_f(k+2,j) + multConst[2] + multConst[1] + auPenalty(k+2, j) + d5 + ps.total() <= mfe_ + delta_) {
                        ps_t ps1(ps);
                        ps1.push(segment(i, k, lM, FM[i][k]));
                        ps1.push(segment(k+2, j, lV, V_f(k+2,j)));
                        ps1.accumulate(multConst[2] + multConst[1] + auPenalty(k+2, j) + d5);
                        ps1.update(k+2,j,'(',')');
                        ps1.update(k+1, d5symb);
                        push_to_gstack(gstack, ps1);
                    }
                    if (FM[i][k] + V_f(k+1,j-1) + multConst[2] + multConst[1] + auPenalty(k+1, j-1) + d3 + ps.total() <= mfe_ + delta_) {
                        ps_t ps1(ps);
                        ps1.push(segment(i, k, lM, FM[i][k]));
                        ps1.push(segment(k+1, j-1, lV, V_f(k+1,j-1)));
                        ps1.accumulate(multConst[2] + multConst[1] + auPenalty(k+1, j-1) + d3);
                        ps1.update(k+1,j-1,'(',')');
                        ps1.update(j, d3symb);
                        push_to_gstack(gstack, ps1);
                    }
                    if (FM[i][k] + V_f(k+2,j-1) + multConst[2] + 2*multConst[1] + auPenalty(k+2, j-1) + d53 + ps.total() <= mfe_ + delta_) {
                        ps_t ps1(ps);
                        ps1.push(segment(i, k, lM, FM[i][k]));
                        ps1.push(segment(k+2, j-1, lV, V_f(k+2,j-1)));
                        ps1.accumulate(multConst[2] + 2*multConst[1] + auPenalty(k+2, j-1) + d53);
                        ps1.update(k+2,j-1,'(',')');
                        ps1.update(k+1, j, d5symb, d3symb);
                        push_to_gstack(gstack, ps1);
                    }
                    break;
                }

            case BOTH_DANGLE:
                {
                    mpq_class d5 = Ed5_new(k+1, j, k);
                    mpq_class d3 = Ed3_new(k+1, j, j+1);
                    mpq_class aup = auPenalty(k+1, j);
                    bonus = d5 + d3 + multConst[2] + aup;
                    if (FM[i][k] + V_f(k+1,j) + bonus + ps.total() <= mfe_ + delta_) {
                        ps_t ps1(ps);
                        ps1.push(segment(i, k, lM, FM[i][k]));
                        ps1.push(segment(k+1, j, lV, V_f(k+1,j)));
                        ps1.accumulate(bonus);
                        ps1.update(k+1,j,'(',')');
                        push_to_gstack(gstack, ps1);
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

            switch (g_dangles) {
            case NO_DANGLE:
                {
                    bonus = multConst[2] + multConst[1]*(k-i+1) + auPenalty(k+1, j);
                    if (V_f(k+1,j) + bonus + ps.total() <= mfe_ + delta_) {
                        ps_t ps1(ps);
                        ps1.push(segment(k+1, j, lV, V_f(k+1,j)));
                        ps1.accumulate(bonus);
                        ps1.update(k+1, j, '(', ')');
                        push_to_gstack(gstack, ps1);
                    }
                    break;
                }

            case CHOOSE_DANGLE:
                {
                    mpq_class d5 = Ed5_pair(k+2, j);
                    mpq_class d3 = Ed3_pair(k+1, j-1);
                    mpq_class d53 = Ed5_pair(k+2, j-1) + Ed3_pair(k+2, j-1);
                    if (V_f(k+1,j) + multConst[2] + multConst[1]*(k+1 - i) + auPenalty(k+1,j) + ps.total() <= mfe_ + delta_) {
                        ps_t ps1(ps);
                        ps1.push(segment(k+1, j, lV, V_f(k+1,j)));
                        ps1.accumulate(multConst[2] + multConst[1]*(k+1 - i) + auPenalty(k+1,j));
                        ps1.update(k+1, j, '(', ')');
                        push_to_gstack(gstack, ps1);
                    }
                    if (V_f(k+2,j) + multConst[2] + multConst[1]*(k+2 - i) + auPenalty(k+2,j) + d5 + ps.total() <= mfe_ + delta_) {
                        ps_t ps1(ps);
                        ps1.push(segment(k+2, j, lV, V_f(k+2,j)));
                        ps1.accumulate(multConst[2] + multConst[1]*(k+2 - i) + auPenalty(k+2,j) + d5);
                        ps1.update(k+2, j, '(', ')');
                        ps1.update(k+1, d5symb);
                        push_to_gstack(gstack, ps1);
                    }
                    if (V_f(k+1,j-1) + multConst[2] + multConst[1]*(k+1 - i + 1) + auPenalty(k+1,j-1) + d3 + ps.total() <= mfe_ + delta_) {
                        ps_t ps1(ps);
                        ps1.push(segment(k+1, j-1, lV, V_f(k+1,j-1)));
                        ps1.accumulate(multConst[2] + multConst[1]*(k+1 - i + 1) + auPenalty(k+1,j-1) + d3);
                        ps1.update(k+1, j-1, '(', ')');
                        ps1.update(j, d3symb);
                        push_to_gstack(gstack, ps1);
                    }
                    if (V_f(k+2,j-1) + multConst[2] + multConst[1]*(k+2 - i + 1) + auPenalty(k+2,j-1) + d53 + ps.total() <= mfe_ + delta_) {
                        ps_t ps1(ps);
                        ps1.push(segment(k+2, j-1, lV, V_f(k+2,j-1)));
                        ps1.accumulate(multConst[2] + multConst[1]*(k+2 - i + 1) + auPenalty(k+2,j-1) + d53);
                        ps1.update(k+2, j-1, '(', ')');
                        ps1.update(k+1, j, d5symb, d3symb);
                        push_to_gstack(gstack, ps1);
                    }
                    break;
                }

            case BOTH_DANGLE:
                {
                    mpq_class d5 = Ed5_new(k+1, j, k);
                    mpq_class d3 = Ed3_new(k+1, j, j+1);
                    bonus = d5 + d3 + multConst[2] + multConst[1]*(k-i+1) + auPenalty(k+1, j);
                    if (V_f(k+1,j) + bonus + ps.total() <= mfe_ + delta_) {
                        ps_t ps1(ps);
                        ps1.push(segment(k+1, j, lV, V_f(k+1,j)));
                        ps1.accumulate(bonus);
                        ps1.update(k+1, j, '(', ')');
                        push_to_gstack(gstack, ps1);
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
    };

    /*
    // Helper method to manage mess in traceWM
    void pushV(ps_t& ps, ps_stack_t& gstack, int i, int j, mpq_class bonus) {
        if (V_f(i, j) + bonus + ps.total() <= mfe_ + delta_) {
            ps_t ps_new(ps);
            ps_new.accumulate(bonus);
            ps_new.push(segment(i, j, lV, V_f(i, j)));
            push_to_gstack(gstack, ps_new);
        }

        if (WMPrime[i][j] + ps.total() <= mfe_ + delta_) {
            ps_t ps_new(ps);
            ps_new.push(segment(i, j, lWMPrime, WMPrime[i][j]));
            push_to_gstack(gstack, ps_new);
        }

        if (WM_f(i+1,j) + multConst[1] + ps.total() <= mfe_ + delta_) {
            ps_t ps_new(ps);
            ps_new.accumulate(multConst[1]);
            ps_new.push(segment(i+1, j, lWM, WM_f(i+1, j)));
            push_to_gstack(gstack, ps_new);
        }

        if (WM_f(i, j-1) + multConst[1] + ps.total() <= mfe_ + delta_) {
            ps_t ps_new(ps);
            ps_new.accumulate(multConst[1]);
            ps_new.push(segment(i, j-1, lWM, WM_f(i, j-1)));
            push_to_gstack(gstack, ps_new);
        }
    }

    void traceWM(int i, int j, ps_t& ps, ps_stack_t& gstack) {
        if (DEBUG) printf("traceWM(%i, %i):\n\tpartial structure %s\n", i, j, ps.str.c_str());

        switch (g_dangles) {
        case NO_DANGLE:
            {
                pushV(ps, gstack, i, j, auPenalty(i, j) + multConst[2]);
                break;
            }

        case CHOOSE_DANGLE:
            {
                mpq_class d5 = Ed5_pair(i, j-1);
                mpq_class d3 = Ed3_pair(i+1, j);
                mpq_class d53 = Ed5_pair(i+1, j-1) + Ed3_pair(i+1, j-1);
                pushV(ps, gstack, i, j, auPenalty(i, j) + multConst[2]);
                pushV(ps, gstack, i+1, j, auPenalty(i+1, j) + multConst[2] + multConst[1] + d5);
                pushV(ps, gstack, i, j-1, auPenalty(i, j-1) + multConst[2] + multConst[1] + d3);
                pushV(ps, gstack, i+1, j-1, auPenalty(i+1, j-1) + multConst[2] + 2*multConst[1] + d53);
                break;
            }

        case BOTH_DANGLE:
            {
                mpq_class d5 = Ed5(j,i,j+1);
                mpq_class d3 = (i==1)?Ed3(j,i,length):Ed3(j,i,i-1);
                pushV(ps, gstack, i, j, auPenalty(i, j) + multConst[2] + d5 + d3);
                break;
            }

        default:
            {
                exit(EXIT_FAILURE);
                break;
            }
        };
    }

    void traceWMPrime(int i, int j, ps_t& ps, ps_stack_t& gstack) {
        if (DEBUG) printf("traceWMPrime(%i, %i):\n\tpartial structure %s\n", i, j, ps.str.c_str());
        for (int h = i+TURN+1 ; h <= j-TURN-2; h++) {
            if (WM_f(i,h-1) + WM_f(h,j) + ps.total() <= mfe_ + delta_) {
                ps_t ps_new(ps);
                ps_new.push(segment(i,h-1, lWM, WM_f(i,h-1)));
                ps_new.push(segment(h,j, lWM, WM_f(h,j)));
                push_to_gstack(gstack, ps_new);
            }
        }
    }

    void traceVM(int i, int j, ps_t& ps, ps_stack_t& gstack) {
        if (DEBUG) printf("traceVM(%i, %i):\n\tpartial structure %s\n", i, j, ps.str.c_str());

        mpq_class bonus = 0;

        switch (g_dangles) {
        case NO_DANGLE:
            {
                bonus = multConst[0] + multConst[2] + auPenalty(i, j);
                if (WMPrime[i+1][j-1] + bonus + ps.total() <= mfe_ + delta_) {
                    ps_t ps_new(ps);
                    ps_new.accumulate(bonus);
                    ps_new.push(segment(i+1,j-1, lWMPrime,WMPrime[i+1][j-1] ));
                    push_to_gstack(gstack, ps_new);
                }
                break;
            }

        case CHOOSE_DANGLE:
            {
                mpq_class d5 = Ed5_pair(i+2, j-1);
                mpq_class d3 = Ed3_pair(i+1, j-2);
                mpq_class d53 = Ed5_pair(i+2, j-2) + Ed3_pair(i+2, j-2);
                if (WMPrime[i+1][j - 1] + multConst[0] + multConst[2] + auPenalty(i,j) + ps.total() <= mfe_ + delta_) {
                    ps_t ps_new(ps);
                    ps_new.accumulate(multConst[0] + multConst[2] + auPenalty(i,j));
                    ps_new.push(segment(i+1,j-1, lWMPrime,WMPrime[i+1][j-1] ));
                    push_to_gstack(gstack, ps_new);
                }
                if (WMPrime[i + 2][j - 1] + multConst[0] + multConst[2] + multConst[1] + auPenalty(i+1,j) + d5 + multConst[1] + ps.total() <= mfe_ + delta_) {
                    ps_t ps_new(ps);
                    ps_new.accumulate(multConst[0] + multConst[2] + multConst[1] + auPenalty(i+1,j) + d5 + multConst[1]);
                    ps_new.push(segment(i+2,j-1, lWMPrime,WMPrime[i+2][j-1] ));
                    ps_new.update(i+1, d5symb);
                    push_to_gstack(gstack, ps_new);
                }
                if (WMPrime[i + 1][j - 2] + multConst[0] + multConst[2] + multConst[1] + auPenalty(i, j-1) + d3 + multConst[1] + ps.total() <= mfe_ + delta_) {
                    ps_t ps_new(ps);
                    ps_new.accumulate(multConst[0] + multConst[2] + multConst[1] + auPenalty(i, j-1) + d3 + multConst[1]);
                    ps_new.push(segment(i+1,j-2, lWMPrime,WMPrime[i+1][j-2]));
                    ps_new.update(j-1, d3symb);
                    push_to_gstack(gstack, ps_new);
                }
                if ( WMPrime[i + 2][j - 2] + multConst[0] + multConst[2] + 2*multConst[1] + auPenalty(i+1, j-1) + d53 + multConst[1]*2 + ps.total() <= mfe_ + delta_) {
                    ps_t ps_new(ps);
                    ps_new.accumulate(multConst[0] + multConst[2] + 2*multConst[1] + auPenalty(i+1, j-1) + d53 + multConst[1]*2);
                    ps_new.push(segment(i+2,j-2, lWMPrime,WMPrime[i+2][j-2] ));
                    ps_new.update(i+1, j-1, d5symb, d3symb);
                    push_to_gstack(gstack, ps_new);
                }
                break;
            }

        case BOTH_DANGLE:
            {
                mpq_class d3 = Ed3(i,j,j-1);
                mpq_class d5 = Ed5(i,j,i+1);
                bonus = multConst[0] + multConst[2] + auPenalty(i, j) + d5 + d3;
                if (WMPrime[i+1][j-1] + bonus + ps.total() <= mfe_ + delta_) {
                    ps_t ps_new(ps);
                    ps_new.accumulate(bonus);
                    ps_new.push(segment(i+1,j-1, lWMPrime,WMPrime[i+1][j-1] ));
                    push_to_gstack(gstack, ps_new);
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
    */

    void push_to_gstack(ps_stack_t& gstack, const ps_t& v) {
        if (DEBUG)
            {
                if (v.empty())
                {
                    printf("Pushing partial structure %s with no segments to stack\n", v.str.c_str());
                }
                else {
                    segment smt = v.top();
                    printf("Pushing partial structure %s with top segment (%i, %i) to stack for trace %s\n", v.str.c_str(), smt.i_, smt.j_, lstr[smt.label_]);
                }
            }
        gflag = 1;
        gstack.push(v);
    }
}
