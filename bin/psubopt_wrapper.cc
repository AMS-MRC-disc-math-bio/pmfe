// Copyright (c) 2014 Andrew Gainer-Dewar.

#include "subopt.h"
#include "pmfe_types.h"
#include <iostream>
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
        ("outdir,o", po::value<std::string>(), "Output structure directory")
        ("delta", po::value<std::string>()->default_value("0"), "Energy delta value")
        ("paramdir,p", po::value<std::string>()->default_value("/usr/local/share/pmfe/Turner99/pmfe"), "Turner99 parameter directory")
        ("multiloop-penalty,a", po::value<std::string>(), "Multiloop penalty parameter")
        ("unpaired-penalty,b", po::value<std::string>(), "Unpaired base penalty parameter")
        ("branch-penalty,c", po::value<std::string>(), "Branching helix penalty parameter")
        ("dummy-scaling,d", po::value<std::string>(), "Dummy scaling parameter")
        ("dangle-model,m", po::value<int>()->default_value(1), "Dangle model")
        ("max-count", po::value<int>()->default_value(-1), "Max number of structures to generate")
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
    fs::path seq_file, output_dir, param_dir;

    seq_file = fs::path(vm["sequence"].as<std::string>());
    param_dir = fs::path(vm["paramdir"].as<std::string>());


    if (vm.count("outdir")) {
        output_dir = fs::path(vm["outfile"].as<std::string>());
    } else {
        output_dir = seq_file;
        output_dir.replace_extension("");
    }

    mpq_class delta = mpq_class(vm["delta"].as<std::string>());

    // Set up the parameter vector
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

    // Set up dangle model
    int dangle_model = vm["dangle-model"].as<int>();

    // Set up max count
    int max_count = vm["max-count"].as<int>();

    pmfe::generate_subopt(seq_file.native(), output_dir.native(), delta, params, param_dir.native(), max_count, dangle_model);
    return(0);
}
