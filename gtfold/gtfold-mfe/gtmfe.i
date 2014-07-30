/* File: gtmfe.i */
%module gtmfe
%include "std_string.i"
%include "include/PolytopeVector.h"

%{
#define SWIG_FILE_WITH_INIT
#include "mfe_main.h"
#include "utils.h"
#include "loader.h"
#include "global.h"
#include "energy.h"
#include "algorithms.h"
#include "constraints.h"
#include "traceback.h"
#include "shapereader.h"
#include "omp.h"
%}

PolytopeVector mfe_main(std::string seq_file, std::string output_file, std::string param_dir, int dangle_model = 1);
