// Copyright (c) 2015 Andrew Gainer-Dewar.

#include <iostream>

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

#include "pmfe_types.h"
#include "subopt.h"

namespace fs = boost::filesystem;
namespace po = boost::program_options;

int main(int argc, char* argv[]) {
    // Set up argument processing
    po::options_description desc("Options");
    desc.add_options()
        ("sequence", po::value<std::string>()->required(), "Sequence file")
        ("outfile,o", po::value<std::string>(), "Output file")
        ("delta", po::value<std::string>()->default_value("0"), "Energy delta value")
        ("multiloop-penalty,a", po::value<std::string>(), "Multiloop penalty parameter")
        ("unpaired-penalty,b", po::value<std::string>(), "Unpaired base penalty parameter")
        ("branch-penalty,c", po::value<std::string>(), "Branching helix penalty parameter")
        ("dummy-scaling,d", po::value<std::string>(), "Dummy scaling parameter")
        ("dangle-model,m", po::value<int>()->default_value(1), "Dangle model")
        ("sorted,s", po::bool_switch(), "Sort results in increasing energy order")
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

    // Process file-related options
    fs::path seq_file(vm["sequence"].as<std::string>());
    fs::path out_file(vm["outfile"].as<std::string>());

    // Set up the parameters
    mpq_class delta = pmfe::get_mpq_from_word(vm["delta"].as<std::string>());

    pmfe::ParameterVector params = pmfe::ParameterVector();

    if (vm.count("multiloop-penalty")) {
        params.multiloop_penalty = pmfe::get_mpq_from_word(vm["multiloop-penalty"].as<std::string>());
    };

    if (vm.count("unpaired-penalty")) {
        params.unpaired_penalty = pmfe::get_mpq_from_word(vm["unpaired-penalty"].as<std::string>());
    };

    if (vm.count("branch-penalty")) {
        params.branch_penalty = pmfe::get_mpq_from_word(vm["branch-penalty"].as<std::string>());
    };

    if (vm.count("dummy-scaling")) {
        params.dummy_scaling = pmfe::get_mpq_from_word(vm["dummy-scaling"].as<std::string>());
    };

    params.canonicalize();

    bool sorted = vm["sorted"].as<bool>();

    // Set up dangle model
    pmfe::dangle_mode dangles = pmfe::convert_to_dangle_mode(vm["dangle-model"].as<int>());

    // Get results
    std::vector<pmfe::RNAStructureWithScore> structures = suboptimal_structures(seq_file, params, dangles, delta, sorted);

    for (int i = 0; i < structures.size(); ++i) {
        std::cout << i << "\t" << structures[i] << std::endl;
    }
}
