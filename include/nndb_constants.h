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
*/

/* Amrita: some of the variables defined in this file are not uesed. */

#ifndef _NNDB_CONSTANTS_H
#define _NNDB_CONSTANTS_H

#include "pmfe_types.h"
#include "thread_pool.h"
#include "rational.h"

#include <boost/filesystem.hpp>
#include <boost/multi_array.hpp>
#include <vector>
#include <map>

namespace pmfe {
    namespace fs = boost::filesystem;
    //TODO: Make abstract
    //TODO: Provide Turner99 instance
    class NNDBConstants {
    public:
        Rational maxpen;
        Rational auend;
        Rational gubonus;
        Rational cint; /* cint, cslope, c3 are used for poly C hairpin loops */
        Rational cslope;
        Rational c3;
        bool gail;
        Rational prelog;
        ParameterVector params;

        std::vector<Rational> poppen;
        std::vector<Rational> multConst; /* for multiloop penalties. */
        std::vector<Rational> inter; /* Contains size penalty for internal loops */
        std::vector<Rational> bulge; /* Contain the size penalty for bulges */
        std::vector<Rational> hairpin; /* Contains the size penalty for hairpin loops */

        std::map<std::string, Rational> tloop;

        //TODO: Rebase anything that should be 1-indexed
        boost::multi_array<Rational, 4> tstkh; /* Terminal mismatch energy used in the calculations of hairpin loops */
        boost::multi_array<Rational, 4> tstki; /* Terminal mismatch energy used in the calculations of internal loops */
        boost::multi_array<Rational, 4> stack; /* Stacking energy used to calculate energy of stack loops */
        boost::multi_array<Rational, 4> dangle; /* Contain dangling energy values */
        boost::multi_array<Rational, 8> iloop22; /* 2*2 internal looops */
        boost::multi_array<Rational, 7> iloop21; /* 2*1 internal loops */
        boost::multi_array<Rational, 6> iloop11; /* 1*1 internal loops */

    NNDBConstants(const ParameterVector params = ParameterVector()):
        params(params),
            poppen(5),
            multConst(3),
            inter(31),
            bulge(31),
            hairpin(31),
            tstkh(boost::extents[4][4][4][4]),
            tstki(boost::extents[4][4][4][4]),
            stack(boost::extents[4][4][4][4]),
            dangle(boost::extents[4][4][4][2]),
            iloop22(boost::extents[5][5][5][5][5][5][5][5]),
            iloop21(boost::extents[5][5][5][5][5][5][5]),
            iloop11(boost::extents[5][5][5][5][5][5])
            {};
    };

    class Turner99: public NNDBConstants {
    public:
        Turner99(SimpleThreadPool& thread_pool, const ParameterVector& params = ParameterVector(), const fs::path& param_dir = "Turner99");

    protected:
        SimpleThreadPool& thread_pool;

        void initMiscValues(const fs::path& param_dir);
        void initLoopValues(const fs::path& param_dir);
        void initTstkhValues(const fs::path& param_dir);
        void initTstkiValues(const fs::path& param_dir);
        void initStackValues(const fs::path& param_dir);
        void initDangleValues(const fs::path& param_dir);
        void initTloopValues(const fs::path& param_dir);
        void initIloop11Values(const fs::path& param_dir);
        void initIloop21Values(const fs::path& param_dir);
        void initIloop22Values(const fs::path& param_dir);
    };
}

#endif
