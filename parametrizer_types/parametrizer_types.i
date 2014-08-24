/* File: parametrizer_types.i */
%module parametrizer_types

%include "parametrizer_types.h"
%include "std_pair.i"
%include "std_vector.i"
%include "std_string.i"

%{
  #define SWIG_FILE_WITH_INIT
  #include "parametrizer_types.h"
%}

%template(pairll) std::pair< long, long >;
%template(vecll)  std::vector< std::pair< long, long > >;
%template(vecstr) std::vector< std::string >;

%extend ScoreVector{
  %pythoncode{
     def get_python_fractions_dict(self):
         from fractions import Fraction
         pairs = self.get_pairs()
         result = {"multiloops" : pairs[0][0],
                   "unpaired": pairs[1][0],
                   "branches": pairs[2][0],
                   "w": Fraction(pairs[3][0], pairs[3][1]),
                   "energy": Fraction(pairs[4][0], pairs[4][1]),
         }
                  
         return result
       }
 }

%extend ParameterVector{
  %pythoncode{
     def get_python_fractions_dict(self):
         from fractions import Fraction
         pairs = self.get_pairs()
         result = {"multiloop penalty" : Fraction(pairs[0][0], pairs[0][1]),
                   "unpaired base penalty": Fraction(pairs[1][0], pairs[1][1]),
                   "branching helix penalty": Fraction(pairs[2][0], pairs[2][1]),
                   "dummy scaling parameter": Fraction(pairs[3][0], pairs[3][1]),
         }
                  
         return result
       }
 }
