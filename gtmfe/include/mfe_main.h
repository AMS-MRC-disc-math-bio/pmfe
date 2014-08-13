#ifndef _MFE_MAIN_H_
#define _MFE_MAIN_H_

#include "helper-structs.h"

PolytopeVector mfe_main(std::string seq_file, std::string output_file, std::string param_dir, long double a=3.4, long double b=0, long double c=0.4, long double d=1, int dangle_model = 1);

void init_fold(const char* seq);
void free_fold(int len);

#endif
