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

#include "mfe.h"
#include "utils.h"
#include "loader.h"
#include "global.h"
#include "energy.h"
#include "algorithms.h"
#include "traceback.h"
#include "parametrizer_types.h"

#include "RNAScoring.h"

#include <gmpxx.h>

using namespace std;

static string seqfile = "";
static string outputFile = "";
static string outputDir = "";

static int dangles=-1;

static int nThreads = -1;

void init_fold(const char* seq, ParameterVector params) {
    assert(seq != NULL);
    int len = strlen(seq);

    init_global_params(len, params);

    if (!encodeSequence(seq)) {
        free_fold(len);
        exit(0);
    }

    create_tables(len);

    g_nthreads = nThreads;
    g_unamode  = false;
    g_mismatch = false;
    g_verbose  = false;
    g_prefilter_mode  = false;
    g_prefilter1  = 2;
    g_prefilter2  = 2;
    g_dangles = dangles;

}

void free_fold(int len) {
    free_tables(len);
    free_global_params();
}

ScoreVector mfe(string seq_file, string output_file, string param_dir, int dangle_model) {
    ParameterVector params = ParameterVector();
    return mfe(seq_file, output_file, param_dir, params, dangle_model);
}

ScoreVector mfe(string seq_file, string output_file, string param_dir, ParameterVector params, int dangle_model) {
    std::string seq;
    mpq_class energy;

    dangles = dangle_model;

    if (!(dangles == 0 || dangles == 1 || dangles == 2)) {
        dangles = -1;
    }

    if (read_sequence_file(seq_file.c_str(), seq) == FAILURE) {
        printf("Failed to open sequence file: %s.\n\n", seqfile.c_str());
        exit(-1);
    }

    init_fold(seq.c_str(), params);

    // Read in thermodynamic parameters.
    readThermodynamicParameters(param_dir.c_str());

    // Compute the minimum free energy
    energy = calculate(seq.length());

    // Find the associated structure
    ScoreVector result;
    result = trace(seq.length());

    if (result.energy != energy) {
        std::cerr << "Energy traceback is inconsistent!" << std::endl;
        std::cerr << params << std::endl;
        std::cerr << result << std::endl;
        std::cerr << "DP energy: " << energy.get_str(10) << std::endl;
        std::cerr << "Traceback energy: " << result.energy.get_str(10) << std::endl << std::endl;
    }

    // Write out the resulting structure
    save_ct_file(output_file, seq, energy);

    // Find the classical energy
    mpq_class classical_energy = rnascoring::get_classical_score(output_file, "rna-scoring/data/Turner99", dangle_model);
    classical_energy.canonicalize();

    // Calculate w
    result.w = classical_energy - (result.multiloops * multiloop_default + result.unpaired * unpaired_default + result.branches * branch_default);
    result.canonicalize();

    // Check that the w calculation produced a consistent result
    mpq_class formula_energy = result.multiloops * params.multiloop_penalty + result.unpaired * params.unpaired_penalty + result.branches * params.branch_penalty + result.w * params.dummy_scaling;
    formula_energy.canonicalize();

    // And alert the user if not
    if (result.energy != formula_energy) {
        std::cerr << "Energy calculation is inconsistent!" << std::endl;
        std::cerr << params << std::endl;
        std::cerr << result << std::endl;
        std::cerr << "Formula energy: " << formula_energy.get_str(10) << std::endl;
        std::cerr << "Classical energy: " << classical_energy.get_str(10) << std::endl << std::endl;
    };

    // Clean up
    free_fold(seq.length());

    return result;
}
