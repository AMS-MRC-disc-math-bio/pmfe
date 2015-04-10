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

#include <gmpxx.h>
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
        int numoftloops;
        mpq_class maxpen;
        mpq_class auend;
        mpq_class gubonus;
        mpq_class cint; /* cint, cslope, c3 are used for poly C hairpin loops */
        mpq_class cslope;
        mpq_class c3;
        bool gail;
        mpq_class prelog;
        mpq_class INFINITY_;

        std::vector<mpq_class> poppen;
        std::vector<mpq_class> eparam;
        std::vector<mpq_class> multConst; /* for multiloop penalties. */
        std::vector<mpq_class> inter; /* Contains size penalty for internal loops */
        std::vector<mpq_class> bulge; /* Contain the size penalty for bulges */
        std::vector<mpq_class> hairpin; /* Contains the size penalty for hairpin loops */

        std::map<std::string, mpq_class> tloop;

        //TODO: Rebase anything that should be 1-indexed
        boost::multi_array<mpq_class, 4> tstkh; /* Terminal mismatch energy used in the calculations of hairpin loops */
        boost::multi_array<mpq_class, 4> tstki; /* Terminal mismatch energy used in the calculations of internal loops */
        boost::multi_array<mpq_class, 4> stack; /* Stacking energy used to calculate energy of stack loops */
        boost::multi_array<mpq_class, 4> dangle; /* Contain dangling energy values */
        boost::multi_array<mpq_class, 8> iloop22; /* 2*2 internal looops */
        boost::multi_array<mpq_class, 7> iloop21; /* 2*1 internal loops */
        boost::multi_array<mpq_class, 6> iloop11; /* 1*1 internal loops */

    NNDBConstants(const ParameterVector params = ParameterVector()):
        poppen(5),
            eparam(11),
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
            {
                INFINITY_ = mpq_class(9999999999999);
                if (abs(params.dummy_scaling > 1)) {
                    INFINITY_ *= abs(params.dummy_scaling);
                }
            };
    };

    class Turner99: public NNDBConstants {
    public:
        Turner99(const ParameterVector params = ParameterVector());

    protected:
        ParameterVector params;
        void initMiscValues();
        void initLoopValues();
        void initTstkValues();
        void initStackValues();
        void initDangleValues();
        void initTloopValues();
        void initIloop11Values();
        void initIloop21Values();
        void initIloop22Values();
    };
}

#endif
