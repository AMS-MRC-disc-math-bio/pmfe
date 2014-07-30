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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cassert>

#include "main.h"
#include "mfe_main.h"
#include "utils.h"
#include "loader.h"
//#include "options.h"
#include "global.h"
#include "energy.h"
#include "algorithms.h"
#include "constraints.h"
#include "traceback.h"
#include "shapereader.h"

using namespace std;

static bool PARAM_DIR = false;
//static bool LIMIT_DISTANCE;
static bool CONS_ENABLED = false;
static bool VERBOSE = false;
static bool SILENT = false;
//static bool SHAPE_ENABLED = false;

//extern int SHAPE_ENABLED;

static bool T_MISMATCH = false;
static bool UNAMODE = false;
static bool RNAMODE = false;
static bool b_prefilter = false;

static string seqfile = "";
static string constraintsFile = "";
static string outputPrefix = "";
static string outputFile = "";
static string energyDecomposeOutFile = "";
static string outputDir = "";
static string shapeFile = "";
static string paramDir; // default value

static int dangles=-1;
static int prefilter1=2;
static int prefilter2=2;
static int print_energy_decompose = 0;

static int nThreads = -1;
static int contactDistance = -1;

void init_fold(const char* seq) {
  assert(seq != NULL);
  int len = strlen(seq);

  init_global_params(len);

  if (!encodeSequence(seq)) {
    free_fold(len);
    exit(0);
  }

  create_tables(len);

  if (CONS_ENABLED) {
    init_constraints(constraintsFile.c_str(), len);
  }

  if (SHAPE_ENABLED) {
    readSHAPEarray(shapeFile.c_str(),len);
  }

  if (UNAMODE) {
    if (T_MISMATCH) if(!SILENT) printf("Ignoring -m option, using --unafold\n");
    if (PARAM_DIR) if(!SILENT) printf("Ignoring -p option, using --unafold\n");
    if (dangles == 0 || dangles == 1 || dangles == 2) 
      if(!SILENT) printf("Ignoring -d option, using --unafold\n");
    if (b_prefilter == 1) 
      if(!SILENT) printf("Ignoring --prefilter option, using --unafold\n");
    T_MISMATCH = false;
    PARAM_DIR = false;
    dangles = -1;
    b_prefilter = false;
  }

  if (RNAMODE) {
    if (T_MISMATCH) if(!SILENT) printf("Ignoring -m option, using --rnafold\n");
    if (PARAM_DIR) if(!SILENT) printf("Ignoring -p option, using --rnafold\n");
    if (dangles == 0 || dangles == 1 || dangles == 2) 
      if(!SILENT) printf("Ignoring -d option, using --rnafold\n");
    if (b_prefilter == 1) 
      if(!SILENT) printf("Ignoring --prefilter option, using --rnafold\n");
    T_MISMATCH = false;
    PARAM_DIR = false;
    dangles = -1;
    b_prefilter = false;
  }

  if ((dangles == 0 || dangles == 1 ||dangles == 2) && !UNAMODE && !RNAMODE) {
    if (T_MISMATCH) if(!SILENT) printf("Ignoring -m option, using -d option\n");
    T_MISMATCH = false;
  } else {
    if (dangles != -1 && !UNAMODE && !RNAMODE) if(!SILENT) printf("Ignoring -d as it accept 0 1 or 2 only\n");	
    dangles = -1;
  }
  if(dangles==1) dangles=-1;

  g_nthreads = nThreads;
  g_unamode  = UNAMODE;
  g_mismatch = T_MISMATCH;
  g_verbose  = VERBOSE;
  g_prefilter_mode  = b_prefilter;
  g_prefilter1  = prefilter1;
  g_prefilter2  = prefilter2;
  g_dangles = dangles;

#ifdef DEBUG
  if(!SILENT) printf("g_nthreads = %d\n", g_nthreads);
  if(!SILENT) printf("g_unamode = %d\n", g_unamode);
  if(!SILENT) printf("g_mismatch = %d\n", g_mismatch);
  if(!SILENT) printf("g_prefilter_mode = %d\n", g_prefilter_mode);
  if(!SILENT) printf("g_dangles = %d\n", g_dangles);

#endif

}

void free_fold(int len) {
	if (CONS_ENABLED) 
		free_constraints(len);
	if (SHAPE_ENABLED){
		free_shapeArray(len);
	}

	free_tables(len);
	free_global_params();
}


int mfe_main(string seq_file, string output_file, string param_dir, int dangle_model) {
        std::string seq;
        int energy;
	
        dangles = dangle_model;
             
        if (!(dangles == 0 || dangles == 1 || dangles == 2)) {
          dangles = -1;
        }

        paramDir = param_dir;
        PARAM_DIR = true;

        outputFile = output_file;

	if (read_sequence_file(seq_file.c_str(), seq) == FAILURE) {
		printf("Failed to open sequence file: %s.\n\n", seqfile.c_str());
		exit(-1);
	}

	init_fold(seq.c_str());
	
	// Read in thermodynamic parameters. Always use Turner99 data (for now)
        readThermodynamicParameters(paramDir.c_str(), PARAM_DIR, UNAMODE, RNAMODE, T_MISMATCH);

	energy = calculate(seq.length()) ;
        printf("%d", energy);
	
	trace(seq.length(), print_energy_decompose, energyDecomposeOutFile.c_str());
	
	save_ct_file(outputFile, seq, energy);

	free_fold(seq.length());
	
        return EXIT_SUCCESS;
}

//double calculate_mfe(int argc, char** argv) {
double calculate_mfe(std::string seq) {
	int energy;
	fflush(stdout);
	double t1 = get_seconds();
	energy = calculate(seq.length()) ; 
	t1 = get_seconds() - t1;
	/*if (energy >= MAXENG)	
	  printf("- Minimum Free Energy: %12.4f kcal/mol\n", 0.00);
	  else
	  printf("- Minimum Free Energy: %12.4f kcal/mol\n", energy/100.00);
	  printf("- MFE runtime: %9.6f seconds\n", t1);*/

	//free_fold(seq.length());
	return energy/100.0;
}
