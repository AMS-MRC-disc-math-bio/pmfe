// Copyright (c) 2014 Andrew Gainer-Dewar.

#include "mfe.h"
#include "pmfe_types.h"

#include <iostream>
#include <omp.h>
#include <string>

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

#define BOOST_LOG_DYN_LINK 1 // Fix an issue with dynamic library loading
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

int main(int argc, char * argv[]) {
    // Set up argument processing
    po::options_description desc("Options");
    desc.add_options()
        ("sequence", po::value<std::string>()->required(), "Sequence file")
        ("verbose,v", po::bool_switch()->default_value(false), "Write verbose debugging output")
        ("multiloop-penalty,a", po::value<std::string>(), "Multiloop penalty parameter")
        ("unpaired-penalty,b", po::value<std::string>(), "Unpaired base penalty parameter")
        ("branch-penalty,c", po::value<std::string>(), "Branching helix penalty parameter")
        ("dummy-scaling,d", po::value<std::string>(), "Dummy scaling parameter")
        ("dangle-model,m", po::value<int>()->default_value(1), "Dangle model")
        ("num-threads,t", po::value<int>()->default_value(0), "Number of threads")
        ("help,h", "Display this help message")
        ;

    po::positional_options_description p;
    p.add("sequence", -1);
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);

    if (vm.count("help") or argc == 1) {
        std::cout << desc << std::endl;
        return 1;
    };

    po::notify(vm);

    // Process thread-related options
    size_t num_threads = (vm["num-threads"].as<int>());
    omp_set_num_threads(num_threads);

    // Process logging-related options
    bool verbose = vm["verbose"].as<bool>();
    if (verbose) {
        boost::log::core::get()->set_filter(
            boost::log::trivial::severity >= boost::log::trivial::debug);
    } else {
        boost::log::core::get()->set_filter
            (boost::log::trivial::severity >= boost::log::trivial::warning);
    }

    // Process file-related options
    fs::path seq_file(vm["sequence"].as<std::string>());

    // Set up the parameter vector
    pmfe::ParameterVector params = pmfe::ParameterVector();

    if (vm.count("multiloop-penalty")) {
        params.multiloop_penalty = pmfe::get_rational_from_word(vm["multiloop-penalty"].as<std::string>());
    };

    if (vm.count("unpaired-penalty")) {
        params.unpaired_penalty = pmfe::get_rational_from_word(vm["unpaired-penalty"].as<std::string>());
    };

    if (vm.count("branch-penalty")) {
        params.branch_penalty = pmfe::get_rational_from_word(vm["branch-penalty"].as<std::string>());
    };

    if (vm.count("dummy-scaling")) {
        params.dummy_scaling = pmfe::get_rational_from_word(vm["dummy-scaling"].as<std::string>());
    };

    params.canonicalize();

    // Setup dangle model
    pmfe::dangle_mode dangles = pmfe::convert_to_dangle_mode(vm["dangle-model"].as<int>());

    pmfe::RNAStructureWithScore result = pmfe::mfe(seq_file, params, dangles);
    std::cout << result << std::endl;
    return(0);
}
