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
#include "constraints.h"
#include "traceback.h"
#include "shapereader.h"
%}

PolytopeVector mfe_main(std::string seq_file, std::string output_file, std::string param_dir, double a=10.1, double b=-.3, double c=-.3, double d=1, int dangle_model = 1);
