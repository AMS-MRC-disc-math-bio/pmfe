// Copyright (c) Andrew Gainer-Dewar 2015

#include <vector>

#include <boost/filesystem.hpp>

#include "pmfe_types.h"
#include "nntm.h"
#include "nndb_constants.h"

namespace fs = boost::filesystem;

namespace pmfe {
    std::vector<RNAStructureWithScore> suboptimal_structures(const fs::path seq_file, const ParameterVector& params, const dangle_mode& dangles, const Rational& delta, bool sorted) {
        Turner99 constants(params);
        RNASequence seq(seq_file);
        NNTM energy_model(constants, dangles);
        RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        std::vector<RNAStructureWithScore> results = energy_model.suboptimal_structures(seq_annotated, delta, sorted);
        return results;
    }
}
