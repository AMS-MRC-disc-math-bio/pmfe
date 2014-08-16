/* File: iB4e.i */
%module(directors="1") iB4e
%feature("director");

%include "euclideanvector.h"
%include "BBpolytope.h"

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
