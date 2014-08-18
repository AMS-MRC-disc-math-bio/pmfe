/* File: iB4e.i */
%module(directors="1") iB4e
%feature("director");

%include "euclideanvector.h"
%include "BBpolytope.h"
%include "std_pair.i"
%include "std_vector.i"

%{
#define SWIG_FILE_WITH_INIT
#include "BBpolytope.h"
#include "config.h"
#include "euclideanvector.h"
#include "faces.h"
#include "linalg.h"
#include "stack.h"
#include "gmpxx.h"
%}

%template(pairll) std::pair<long, long>;
%template(vecEV) std::vector<EuclideanVector>;

// This is easier than writing a TypeMap
%extend EuclideanVector{
  %pythoncode{
    def get_mpq_values(self):
        import gmpy2
        split_values = self.get_split_values()
        return [gmpy2.mpq(val[0], val[1]) for val in split_values]

    def get_split_values(self):
        split_values = [self.get_split_value(i) for i in xrange(self.dimension)]
        return split_values

    def set_mpq_values(self, mpq_values):
        import gmpy2
        split_values = [(long(value.numerator), long(value.denominator)) for value in mpq_values]
        self.set_split_values(split_values)                      

    def set_split_values(self, split_values):
        for index in xrange(self.dimension):
            val = split_values[index]
            self.set_split_value(index, val[0], val[1])
   }
 }


%extend BBPolytope{
  %pythoncode{
    def get_split_vertices(self):
        return [vertex.get_split_values() for vertex in self.vertices]

    def get_mpq_vertices(self):
        return [vertex.get_mpq_values() for vertex in self.vertices]
  }
 }
