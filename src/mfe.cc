/*
  GTfold: compute minimum free energy of RNA secondary structure
  Copyright (C) 2007-2011  David A. Bader, Christine E. Heitsch,
  and Steve C. Harvey
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

#include <fstream>
#include <string>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cassert>
#include <iostream>

#include <gmpxx.h>
#include <boost/filesystem.hpp>

#include "mfe.h"
#include "nndb_constants.h"
#include "nntm.h"
#include "pmfe_types.h"
#include "thread_pool.h"

namespace pmfe {
    namespace fs = boost::filesystem;

    ScoreVector mfe_pywrap(std::string seq_file, ParameterVector params, int dangle_model) {
        return mfe(seq_file, params, convert_to_dangle_mode(dangle_model));
    }

    ScoreVector mfe(fs::path seq_file, fs::path param_dir, dangle_mode dangles) {
        ParameterVector params = ParameterVector();
        return mfe(seq_file, params, dangles);
    }

    ScoreVector mfe(fs::path seq_file, ParameterVector params, dangle_mode dangles) {
        // Construct a thread pool
        SimpleThreadPool thread_pool;

        std::cout << "Reading constants...";
        // Read in thermodynamic parameters.
        Turner99 constants(params);
        std::cout << "done." << std::endl;

        std::cout << "Reading sequence...";
        // Read in the sequence
        RNASequence seq (seq_file);
        std::cout << seq << " done." << std::endl;

        // Compute the minimum free energy
        std::cout << "Setting up energy model...";
        NNTM energy_model(constants, dangles, thread_pool);
        std::cout << "done." << std::endl;

        std::cout << "Computing energy tables...";
        RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        std::cout << "done." << std::endl;

        std::cout << "Computing MFE...";
        mpq_class energy = energy_model.minimum_energy(seq_annotated);
        std::cout << "done. Energy: " << energy.get_str() << std::endl;

        // Find the associated structure
        std::cout << "Computing MFE structure...";
        RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);
        std::cout << "done. Structure: " << scored_structure << std::endl;

        return scored_structure.score;
    }
}
