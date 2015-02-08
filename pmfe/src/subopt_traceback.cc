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
#include <gmpxx.h>

namespace pmfe {
    using std::pair;
    using std::cout;
    using std::endl;

    const char* lstr[] = {"W", "V", "VBI", "VM", "WM", "WMPrime", "fm", "fm1"};

    void (*trace_func[8]) (int i, int j, ps_t& ps, ps_stack_t& gs);
    static mpq_class mfe_ = INFINITY_;
    static mpq_class delta_ = 0;
    static int length = -1;
    static int gflag = 0;

    static mpq_class FM1[1500][1500] = {{mpq_class(0)}};
    static mpq_class FM[1500][1500] = {{mpq_class(0)}};

    static inline mpq_class Ed5_new(int i, int j, int k) {
        return (k!=0) ? Ed3(j, i, k) : Ed3(j,i,length);
    }

    static inline mpq_class Ed3_new(int i, int j, int k) {
        return (k!=length+1)?Ed5(j, i, k) : Ed5(j, i, 1);
    }

    void calculate_fm1() {

        for (int i = 1; i <= length; ++i) {
            for (int j = i+1; j <= length; ++j) {
                mpq_class min = INFINITY_;
                for (int l = i+TURN+1; l <= j; ++l) {
                    mpq_class d5 = Ed5_new(i,l,i-1);
                    mpq_class d3 = Ed3_new(i,l,l+1);

                    mpq_class fm1 = 0, bonus = 0;
                    switch (g_dangles) {
                    case NO_DANGLE:
                        bonus = auPenalty(i,l) + multConst[1]*(j-l) + multConst[2];
                        break;

                    case BOTH_DANGLE:
                        bonus = auPenalty(i,l) + d5 + d3 + multConst[1]*(j-l) + multConst[2];
                        break;

                    default:
                        exit(EXIT_FAILURE);
                        break;
                    }
                    fm1 = V_f(i,l) + bonus;
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
        char buff[4096];

        int count = 0;
        length = len;

        calculate_fm1();
        calculate_fm();

        ps_stack_t gstack;

        // initialize the partial structure, segment stack = {[1,n]}, label = W, list_bp = {}
        ps_t first(0, len);
        first.push(segment(1, len, lW, W[len]));
        gstack.push(first); // initialize the partial structure stack

        while (1) {
            if (gstack.empty()) break; // exit
            ps_t ps = gstack.top();
            gstack.pop();

            if (ps.empty()) {
                count++;
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
                    (*trace_func[smt.label_])(smt.i_, smt.j_, ps, gstack);
                }

                // discarded current segment, using remaining ones
                if (!gflag) {
                    ps_t ps1(ps);
                    gstack.push(ps1);
                }
            }
        }
        outfile.close();
        printf("Counts of structure generated=%d\n", count);

#ifdef DEBUG
        //printf("# SS = %d\n", count);
#endif
    }

    ss_map_t subopt_traceback(int len, mpq_class delta, string suboptFile,  int max_structure_count) {
        trace_func[lW] = traceW;
        trace_func[lV] = traceV;
        trace_func[lVBI] = traceVBI;
        trace_func[lVM] = traceVM;
        trace_func[lWM] = traceWM;
        trace_func[lWMPrime] = traceWMPrime;
        trace_func[lM] = traceM;
        trace_func[lM1] = traceM1;

        mfe_ = W[len];
        delta_ = delta;

        ss_map_t subopt_data;
        process(subopt_data, len, suboptFile, max_structure_count);

        return subopt_data;
    }

    void traceV(int i, int j, ps_t& ps, ps_stack_t& gstack) {

        // Hairpin Loop
        if (eH(i,j) + ps.total()  <= mfe_ + delta_) {
            ps_t ps1(ps);
            ps1.accumulate(eH(i,j));
            ps1.update(i, j, '(', ')');
            push_to_gstack(gstack, ps1);
        }

        // Stack
        if (eS(i, j) + V_f(i+1, j-1) + ps.total() <= mfe_ + delta_) {
            ps_t ps1(ps);
            ps1.push(segment(i+1, j-1, lV, V_f(i+1, j-1)));
            ps1.accumulate(eS(i,j));
            ps1.update(i, j , '(', ')');
            push_to_gstack(gstack, ps1);
        }

        // Internal Loop
        if (VBI_f(i,j) + ps.total() <= mfe_ + delta_) {
            traceVBI(i,j,ps,gstack);
        }

        // Multiloop

        // UNIQUE version
        int k;

        for (k = i+2; k <= j-TURN-1; ++k) {

            mpq_class kenergy1 = FM[i+1][k] + FM1[k+1][j-1];
            mpq_class d5 = Ed5(i, j, i+1);
            mpq_class d3 = Ed3(i, j, j-1);
            mpq_class aup = auPenalty(i,j);

            mpq_class kenergy2 = 0;

            switch (g_dangles) {
            case NO_DANGLE:
                kenergy2 = aup + multConst[0] + multConst[2];
                break;

            case BOTH_DANGLE:
                kenergy2 = d5 + d3 + aup + multConst[0] + multConst[2];
                break;

            default:
                        exit(EXIT_FAILURE);
                        break;
            };

            mpq_class kenergy_total = kenergy1 + kenergy2;
            if (kenergy_total + ps.total() <= mfe_ + delta_) {
                ps_t ps1(ps);
                ps1.push(segment(i+1,k, lM, FM[i+1][k]));
                ps1.push(segment(k+1,j-1, lM1, FM1[k+1][j-1]));
                ps1.accumulate(kenergy2);
                ps1.update(i,j,'(',')');
                push_to_gstack(gstack, ps1);
            }
        }

        // Non-UNIQUE version
        /*
        if (VM_f(i,j) + ps.total() <= mfe_ + delta_) {
            ps_t ps1(ps);
            ps1.push(segment(i, j, lVM, VM_f(i,j)));
            ps1.update(i, j, '(', ')');
            push_to_gstack(gstack, ps1);
        }
        */

    }

    void traceVBI(int i, int j, ps_t& ps, ps_stack_t& gstack) {
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

        for (int l = i; l < j-TURN; ++l) {
            mpq_class wim1 =  std::min(W[l-1], mpq_class(0));
            mpq_class d3 = (l>i)?Ed3(j,l,l-1):mpq_class(0);
            mpq_class d5 = (j<length)?Ed5(j,l,j+1):mpq_class(0);

            mpq_class bonus = 0;
            switch (g_dangles){
            case NO_DANGLE:
                bonus = auPenalty(l, j);
                break;

            case BOTH_DANGLE:
                bonus = auPenalty(l, j) + d3 + d5;
                break;

            default:
                exit(EXIT_FAILURE);
                break;
            };

            if (V_f(l,j) + bonus + wim1 + ps.total() <= mfe_ + delta_ ) {
                ps_t ps1(ps);
                ps1.push(segment(l, j, lV, V_f(l,j)));
                if (wim1 == W[l-1]) ps1.push(segment(i, l-1, lW, W[l-1]));
                ps1.accumulate(bonus);
                push_to_gstack(gstack, ps1);
            }
        }

        if (W[j-1] + ps.total() <= mfe_ + delta_) {
            ps_t ps1(ps);
            ps1.push(segment(i, j-1, lW, W[j-1]));
            push_to_gstack(gstack, ps1);
        }
    }

    void traceM1(int i, int j, ps_t& ps, ps_stack_t& gstack) {

        if (FM1[i][j-1] + multConst[1] + ps.total() <= mfe_ + delta_) {
            ps_t ps1(ps);
            ps1.push(segment(i, j-1, lM1, FM1[i][j-1]));
            ps1.accumulate(multConst[1]);
            push_to_gstack(gstack, ps1);
        }

        mpq_class d5 = Ed5_new(i, j, i-1);
        mpq_class d3 = Ed3_new(i, j, j+1);
        mpq_class aup = auPenalty(i,j);

        mpq_class bonus = 0;

        switch (g_dangles) {
        case NO_DANGLE:
            bonus = aup + multConst[2];
            break;

        case BOTH_DANGLE:
            bonus = d5 + d3 + aup + multConst[2];
            break;

        default:
            exit(EXIT_FAILURE);
            break;
        };

        if ( V_f(i,j) + bonus + ps.total() <= mfe_ + delta_) {
            ps_t ps1(ps);
            ps1.push(segment(i, j, lV, V_f(i,j)));
            ps1.accumulate(bonus);
            ps1.update(i,j,'(',')');
            push_to_gstack(gstack, ps1);
        }
    }

    void traceM(int i, int j, ps_t& ps, ps_stack_t& gstack) {

        if (FM[i][j-1] + multConst[1] + ps.total() <= mfe_ + delta_) {
            ps_t ps1(ps);
            ps1.push(segment(i, j-1, lM, FM[i][j-1]));
            ps1.accumulate(multConst[1]);
            push_to_gstack(gstack, ps1);
        }

        mpq_class d5 = Ed5_new(i, j, i-1);
        mpq_class d3 = Ed3_new(i, j, j+1);
        mpq_class aup = auPenalty(i,j);

        mpq_class bonus = 0;

        switch (g_dangles) {
        case NO_DANGLE:
            bonus = multConst[2] + aup;
            break;

        case BOTH_DANGLE:
            bonus = d5 + d3 + multConst[2] + aup;
            break;

        default:
            exit(EXIT_FAILURE);
            break;
        };

        if (V_f(i,j) + bonus + ps.total() <= mfe_ + delta_) {
            ps_t ps1(ps);
            ps1.push(segment(i, j, lV, V_f(i,j)));
            ps1.accumulate(bonus);
            ps1.update(i,j,'(',')');
            push_to_gstack(gstack, ps1);
        }


        for (int k = i+TURN+1; k <= j-TURN-1; ++k) {

            d5 = Ed5_new(k+1, j, k);
            d3 = Ed3_new(k+1, j, j+1);
            aup = auPenalty(k+1, j);

            mpq_class bonus = 0;

            switch (g_dangles) {
            case NO_DANGLE:
                bonus = multConst[2] + aup;
                break;

            case BOTH_DANGLE:
                bonus = d5 + d3 + multConst[2] + aup;
                break;

            default:
                exit(EXIT_FAILURE);
                break;
            };

            if (FM[i][k] + V_f(k+1,j) + bonus + ps.total() <= mfe_ + delta_) {
                ps_t ps1(ps);
                ps1.push(segment(i, k, lM, FM[i][k]));
                ps1.push(segment(k+1, j, lV, V_f(k+1,j)));
                ps1.accumulate(bonus);
                ps1.update(k+1,j,'(',')');
                push_to_gstack(gstack, ps1);
            }
        }

        for (int k = i; k <= j-TURN-1; ++k) {

            d5 = Ed5_new(k+1, j, k);
            d3 = Ed3_new(k+1, j, j+1);
            aup = auPenalty(k+1, j);

            mpq_class bonus = 0;

            switch (g_dangles) {
            case NO_DANGLE:
                bonus = multConst[2] + multConst[1]*(k-i+1) + aup;
                break;

            case BOTH_DANGLE:
                bonus = d5 + d3 + multConst[2] + multConst[1]*(k-i+1) + aup;
                break;

            default:
                exit(EXIT_FAILURE);
                break;
            };

            if (V_f(k+1,j) + bonus + ps.total() <= mfe_+delta_) {
                ps_t ps1(ps);
                ps1.push(segment(k+1, j, lV, V_f(k+1,j)));
                ps1.accumulate(bonus);
                ps1.update(k+1, j, '(', ')');
                push_to_gstack(gstack, ps1);
            }
        }
    }

    void traceWM(int i, int j, ps_t& ps, ps_stack_t& gstack) {
        mpq_class d3 = (i==1)?Ed3(j,i,length):Ed3(j,i,i-1);
        mpq_class d5 = Ed5(j,i,j+1);

        mpq_class bonus = 0;

        switch (g_dangles) {
        case NO_DANGLE:
            bonus = auPenalty(i, j) + multConst[2];
            break;

        case BOTH_DANGLE:
            bonus = auPenalty(i, j) + multConst[2] + d3 + d5;
            break;

        default:
            exit(EXIT_FAILURE);
            break;
        };

        if (V_f(i,j) + bonus + ps.total() <= mfe_ + delta_) {
            ps_t ps_new(ps);
            ps_new.accumulate(bonus);
            ps_new.push(segment(i,j, lV, V_f(i,j)));
            push_to_gstack(gstack, ps_new);
        }

        if (WMPrime[i][j] + ps.total() <= mfe_ + delta_) {
            ps_t ps_new(ps);
            ps_new.push(segment(i,j, lWMPrime, WMPrime[i][j]));
            push_to_gstack(gstack, ps_new);
        }

        if (WM_f(i+1,j) + multConst[1] + ps.total() <= mfe_ + delta_) {
            ps_t ps_new(ps);
            ps_new.accumulate(multConst[1]);
            ps_new.push(segment(i+1,j, lWM, WM_f(i+1,j)));
            push_to_gstack(gstack, ps_new);
        }

        if (WM_f(i,j-1) + multConst[1] + ps.total() <= mfe_ + delta_) {
            ps_t ps_new(ps);
            ps_new.accumulate(multConst[1]);
            ps_new.push(segment(i,j-1, lWM, WM_f(i,j-1)));
            push_to_gstack(gstack, ps_new);
        }
    }

    void traceWMPrime(int i, int j, ps_t& ps, ps_stack_t& gstack) {
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

        mpq_class d3 = Ed3(i,j,j-1);
        mpq_class d5 = Ed5(i,j,i+1);

        mpq_class bonus = 0;

        switch (g_dangles) {
        case NO_DANGLE:
            bonus = multConst[0] + multConst[2] + auPenalty(i, j);
            break;

        case BOTH_DANGLE:
            bonus = multConst[0] + multConst[2] + auPenalty(i, j) + d3 + d5;
            break;

        default:
            exit(EXIT_FAILURE);
            break;
        };

        if (WMPrime[i+1][j-1] + bonus + ps.total() <= mfe_ + delta_) {
            ps_t ps_new(ps);
            ps_new.accumulate(bonus);
            ps_new.push(segment(i+1,j-1, lWMPrime,WMPrime[i+1][j-1] ));
            push_to_gstack(gstack, ps_new);
        }
    }

    void push_to_gstack(ps_stack_t& gstack, const ps_t& v) {
        gflag = 1;
        gstack.push(v);
    }
}
