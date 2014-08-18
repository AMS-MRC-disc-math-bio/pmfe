#ifndef _MFE_MAIN_H_
#define _MFE_MAIN_H_

#include "helper-structs.h"
#include <gmpxx.h>

PolytopeVector mfe_main(std::string seq_file, std::string output_file, std::string param_dir, mpq_class a=mpq_class(34,10), mpq_class b=mpq_class(0,1), mpq_class c=mpq_class(4,10), mpq_class d=mpq_class(1,1), int dangle_model = 1);

void init_fold(const char* seq);
void free_fold(int len);

#endif
