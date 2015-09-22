#include <vector>
#include "pmfe_types.h"

#ifndef _SUBOPT_H_
#define _SUBOPT_H_

namespace pmfe{
    std::vector<RNAStructureWithScore> suboptimal_structures(const fs::path seq_file, const ParameterVector& params, const dangle_mode& dangles, const Rational& delta, bool sorted = false);
}
#endif
