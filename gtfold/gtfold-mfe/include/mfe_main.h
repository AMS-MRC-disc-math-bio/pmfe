#ifndef _MFE_MAIN_H_
#define _MFE_MAIN_H_

struct PolytopeVector {
  double x, y, z, w, energy;
};

PolytopeVector mfe_main(std::string seq_file, std::string output_file, std::string param_dir, int dangle_model = 1);
void init_fold(const char* seq);
void free_fold(int len);

#endif
