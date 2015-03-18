/*
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

#ifndef _SUBOPT_TRACEBACK_H
#define _SUBOPT_TRACEBACK_H

#include <cassert>
#include <iostream>
#include <sstream>
#include <stack>
#include <map>
#include <vector>
#include <string>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <gmpxx.h>
#include <assert.h>

namespace pmfe{

    enum label {lW=0, lV, lVBI, lM, lM1};
    extern const char* lstr[];

    const char d3symb = '>';
    const char d5symb = '<';
    const char p3symb = ')';
    const char p5symb = '(';
    const char blanksymb = '.';

    struct segment
    {
        int i_;
        int j_;
        int label_;
        mpq_class en_;

        segment(int i, int j, enum label lbl, mpq_class en)
        {
            i_= i; j_ = j;
            label_ = lbl;
            en_= en;
        }

        segment(const segment& ss)
        {
            i_ = ss.i_;
            j_ = ss.j_;
            label_ = ss.label_;
            en_ = ss.en_;
        }

        segment& operator = (const segment& ss)
        {
            if (this != &ss)
            {
                i_ = ss.i_;
                j_ = ss.j_;
                label_ = ss.label_;
                en_ = ss.en_;
            }
            return *this;
        }

            friend std::ostream& operator<<(std::ostream& out, const segment& seg)
        {
            out << "<(" << seg.i_ << ',' << seg.j_ << ')' << " "
                << lstr[seg.label_] << " dG=" << seg.en_ << ">";
            return out;
        }
    };

    typedef segment SEG;
    typedef std::stack<segment>  SEGSTACK;

    struct pstruct
    {
        std::string str;
        SEGSTACK st_segment;
        SEGSTACK st_v; /* used for backtracking in traceWM */

        mpq_class ae_;
        mpq_class le_;

        mpq_class total() const { return ae_ + le_; }

        pstruct() {}

        pstruct(const pstruct& ps)
        {
            str = ps.str;
            st_segment = ps.st_segment;
            st_v = ps.st_v;
            ae_ = ps.ae_;
            le_ = ps.le_;
        }

        pstruct& operator = (const pstruct& ps)
        {
            if (&ps != this)
            {
                str = ps.str;
                st_segment = ps.st_segment;
                st_v = ps.st_v;
                ae_ = ps.ae_;
                le_ = ps.le_;
            }
            return *this;
        }


        void update(int i, int j, char c1, char c2)
            {
                str[i-1] = c1; str[j-1] = c2;
            }

        void update(int i, char c)
            {
                str[i-1] = c;
            }

        void mark_pair(int i, int j)
            {
                assert (str[i-1] == blanksymb);
                str[i-1] = p5symb;

                assert (str[j-1] == blanksymb);
                str[j-1] = p3symb;
            }

        void mark_d5(int i)
            {
                assert (str[i-1] == blanksymb);
                str[i-1] = d5symb;
            }

        void mark_d3(int i)
            {
                assert (str[i-1] == blanksymb);
                str[i-1] = d3symb;
            }

    pstruct(mpq_class ae, int len) : ae_(ae)
        {
            le_ = 0;
            str = std::string(len, blanksymb);
        }

        void accumulate(mpq_class en)
        {
            ae_ += en;
        }

        void push(segment seg)
        {
            st_segment.push(seg);
            le_ += seg.en_;
        }

        void push_v(segment seg)
        {
            st_v.push(seg);
            le_ += seg.en_;
        }

        void pop()
        {
            segment seg = st_segment.top();
            le_ -= seg.en_;
            st_segment.pop();
        }

        void pop_v()
        {
            segment seg = st_v.top();
            le_ -= seg.en_;
            st_v.pop();
        }

        segment top() const
        {
            return st_segment.top();
        }


        segment top_v() const
        {
            return st_v.top();
        }

        bool empty() const
        {
            return st_segment.empty();
        }

        bool empty_v() const
        {
            return st_v.empty();
        }

        void print() const
        {
            SEGSTACK st = st_segment;
            std::cout <<'[' << ' ' ;
            while (!st.empty())
            {
                segment seg = st.top();
                std::cout << seg << ' ';
                st.pop();
            }

            st = st_v;
            std::cout << ']' << ' ' << '[' ;
            while (!st.empty())
            {
                segment seg = st.top();
                std::cout << seg << ' ';
                st.pop();
            }
            std::cout << ']' << ' ' << str << ' ' ;
            std::cout << " ae=" << ae_ << " le=" << le_ << " te=" << ae_+le_  ;
        }
    };

    /* partial structure */
    typedef pstruct ps_t;

    /* partial structure stack */
    typedef std::stack<pstruct> ps_stack_t;

    /*  partial structure map */
    typedef std::map<std::string, pstruct> ps_map_t;

    /* secondary structure map */
    typedef std::map<std::string, mpq_class> ss_map_t;

    void push_to_gstack(ps_stack_t & gs, const ps_t& v);

    ss_map_t subopt_traceback(int len, mpq_class delta, std::string suboptFile, int max_structure_count);

    void subopt_VMhelper(int i, int k, int j, mpq_class dangle_bonus, ps_t& ps, ps_stack_t& gstack);

    void subopt_traceV(int i, int j, ps_t & ps, ps_stack_t & gs);
    void subopt_traceVBI(int i, int j, ps_t & ps, ps_stack_t & gs);
    void subopt_traceW(int i, int j, ps_t & ps, ps_stack_t & gs);

    void subopt_traceM(int i, int j, ps_t & ps, ps_stack_t & gs);
    void subopt_traceM1(int i, int j, ps_t & ps, ps_stack_t & gs);
}
#endif
