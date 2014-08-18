/* File: gtmfe.i */
%module gtmfe
%include "std_string.i"
%include "include/helper-structs.h"

%{
#define SWIG_FILE_WITH_INIT
#include "mfe_main.h"
#include "utils.h"
#include "loader.h"
#include "global.h"
#include "energy.h"
#include "algorithms.h"
#include "traceback.h"
#include "constraints.h"
#include "omp.h"
%}

PolytopeVector mfe_main(std::string seq_file, std::string output_file, std::string param_dir, double a=3.4, double b=0.0, double c=0.4, double d=1, int dangle_model = 1);
