#ifndef _MFE_MAIN_H_
#define _MFE_MAIN_H_

int mfe_main(std::string seq_file, std::string output_file, std::string param_dir, int dangle_model = 1);
//double calculate_mfe(int argc, char** argv);
double calculate_mfe(std::string seq);
void init_fold(const char* seq);
void free_fold(int len);

#endif
