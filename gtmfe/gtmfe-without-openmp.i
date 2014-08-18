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
#include <gmpxx.h>
%}

PolytopeVector mfe_main(std::string seq_file, std::string output_file, std::string param_dir, mpq_class a=3.4, mpq_class b=0.0, mpq_class c=0.4, mpq_class d=1, int dangle_model = 1);
