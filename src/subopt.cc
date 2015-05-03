// Copyright (c) Andrew Gainer-Dewar 2015

#include <vector>

#include <boost/filesystem.hpp>

#include "pmfe_types.h"
#include "nntm.h"
#include "nndb_constants.h"
#include "thread_pool.h"

namespace fs = boost::filesystem;

namespace pmfe {
    std::vector<RNAStructureWithScore> suboptimal_structures(const fs::path seq_file, const ParameterVector& params, const dangle_mode& dangles, const mpq_class& delta, bool sorted) {
        SimpleThreadPool thread_pool;

        Turner99 constants(params);
        RNASequence seq(seq_file);
        NNTM energy_model(constants, dangles, thread_pool);
        RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        std::vector<RNAStructureWithScore> results = energy_model.suboptimal_structures(seq_annotated, delta, sorted);
        return results;
    }
}
