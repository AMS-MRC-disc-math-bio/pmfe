/* File: gtmfe.i */
%module gtmfe
%include "std_string.i"
%include "std_pair.i"
%include "std_vector.i"
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

%template(pairll) std::pair<long, long>;
%template(vecll) std::vector< std::pair<long, long> >;

PolytopeVector mfe_main(std::string seq_file, std::string output_file, std::string param_dir, ParameterVector params, int dangle_model = 1);

PolytopeVector mfe_main(std::string seq_file, std::string output_file, std::string param_dir, int dangle_model = 1);

%extend PolytopeVector{    
  %pythoncode{
    def get_python_numbers(self):
    from fractions import Fraction
      pairs = self.get_pairs()
      result = {"multiloops" : pairs[0][0],
                "branches": pairs[1][0],
                "unpaired": pairs[2][0],
                "w": Fraction(pairs[3][0], pairs[3][1]),
                "energy": Fraction(pairs[4][0], pairs[4][1]),
    }
                  
    return result
      }
}
