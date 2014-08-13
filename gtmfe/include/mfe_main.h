#ifndef _MFE_MAIN_H_
#define _MFE_MAIN_H_

#include "helper-structs.h"

PolytopeVector mfe_main(std::string seq_file, std::string output_file, std::string param_dir, float a=10.1, float b=-.3, float c=-.3, float d=1, int dangle_model = 1);
void init_fold(const char* seq);
void free_fold(int len);

#endif
