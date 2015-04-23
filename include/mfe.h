#ifndef _MFE_MAIN_H_
#define _MFE_MAIN_H_

#include "pmfe_types.h"
#include <gmpxx.h>

#include <boost/filesystem.hpp>

namespace pmfe{
    namespace fs = boost::filesystem;

    ScoreVector mfe(fs::path seq_file, ParameterVector params, dangle_mode dangles = BOTH_DANGLE);
    ScoreVector mfe(fs::path seq_file, dangle_mode dangles = BOTH_DANGLE);

    ScoreVector mfe_pywrap(std::string seq_file, ParameterVector params, int dangle_model = 1);

    void init_fold(const char* seq, ParameterVector params, dangle_mode dangles);
    void free_fold(int len);
}

#endif
