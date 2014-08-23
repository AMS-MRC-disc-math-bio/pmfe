#ifndef _MFE_MAIN_H_
#define _MFE_MAIN_H_

#include "parametrizer-types.h"
#include <gmpxx.h>

ScoreVector mfe_main(std::string seq_file, std::string output_file, std::string param_dir, ParameterVector params, int dangle_model = 1);

ScoreVector mfe_main(std::string seq_file, std::string output_file, std::string param_dir, int dangle_model = 1);

void init_fold(const char* seq);
void free_fold(int len);

#endif
