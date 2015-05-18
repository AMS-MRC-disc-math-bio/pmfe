// Copyright (c) 2014 Andrew Gainer-Dewar.

#include "BBPolytope.h"
#include "nndb_constants.h"
#include "nntm.h"
#include "pmfe_types.h"
#include "rna_polytope.h"
#include "thread_pool.h"

#include <string>

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

int main(int argc, char * argv[]) {
    // Set up argument processing
    po::options_description desc("Options");
    desc.add_options()
        ("sequence", po::value<std::string>()->required(), "Sequence file")
        ("outfile,o", po::value<std::string>(), "Output file")
        ("dangle-model,m", po::value<int>()->default_value(1), "Dangle model")
        ("num-threads,t", po::value<int>()->default_value(0), "Number of threads")
        ("help,h", "Display this help message")
        ;

    po::positional_options_description p;
    p.add("sequence", -1);
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);

    if (vm.count("help") or argc == 1 ) {
        std::cout << desc << std::endl;
        return 1;
    };

    po::notify(vm);

    // Set up thread pool
    size_t num_threads = (vm["num-threads"].as<int>());
    pmfe::SimpleThreadPool thread_pool(num_threads);

    // Set up dangle model
    pmfe::dangle_mode dangles = pmfe::convert_to_dangle_mode(vm["dangle-model"].as<int>());

    // Set up sequence
    fs::path seq_file (vm["sequence"].as<std::string>());
    pmfe::RNASequence sequence(seq_file);

    pmfe::RNAPolytope poly(sequence, dangles, thread_pool);
    poly.build();

    poly.print_statistics();

    fs::path poly_file;
    if (vm.count("outfile")) {
        poly_file = fs::path(vm["outfile"].as<std::string>());
    } else {
        poly_file = seq_file;
        poly_file.replace_extension(".rnapoly");
    }

    poly.write_to_file(poly_file);

    return 0;
}
