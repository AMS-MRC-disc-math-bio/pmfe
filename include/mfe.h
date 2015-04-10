#ifndef _MFE_MAIN_H_
#define _MFE_MAIN_H_

#include "pmfe_types.h"
#include <gmpxx.h>

#include <boost/filesystem.hpp>

namespace pmfe{
    namespace fs = boost::filesystem;

    ScoreVector mfe(fs::path seq_file, fs::path output_file, ParameterVector params, fs::path param_dir = "/usr/local/share/pmfe/Turner99", dangle_mode dangles = BOTH_DANGLE);
    ScoreVector mfe(fs::path seq_file, fs::path output_file, fs::path param_dir = "/usr/local/share/pmfe/Turner99", dangle_mode dangles = BOTH_DANGLE);

    ScoreVector mfe_pywrap(std::string seq_file, std::string output_file, ParameterVector params, std::string param_dir = "/usr/local/share/pmfe/Turner99", int dangle_model = 1);

    void init_fold(const char* seq, ParameterVector params, dangle_mode dangles);
    void free_fold(int len);
}

#endif
