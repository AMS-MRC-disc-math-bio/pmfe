/*
  File: gtmfe-shared.i
  Shared interfaces for gtmfe-SWIG
*/
%module gtmfe
%include "std_string.i"
%include "std_pair.i"
%include "std_vector.i"
%include "../parametrizer-types/parametrizer-types.i"

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
#include "parametrizer-types.h"
#include <gmpxx.h>
%}

ScoreVector mfe_main(std::string seq_file, std::string output_file, std::string param_dir, ParameterVector params, int dangle_model = 1);

ScoreVector mfe_main(std::string seq_file, std::string output_file, std::string param_dir, int dangle_model = 1);
