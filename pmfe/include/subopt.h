#include <string>
#include "pmfe_types.h"

#ifndef _SUBOPT_H_
#define _SUBOPT_H_

namespace pmfe{
    void generate_subopt(std::string seq_file, std::string output_dir, mpq_class delta = 0, ParameterVector params = ParameterVector(), std::string param_dir = "/usr/local/share/pmfe/Turner99/pmfe", int max_structure_count = -1, int dangle_model = 1);
        }
#endif
