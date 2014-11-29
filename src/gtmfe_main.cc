// Copyright (c) 2014 Andrew Gainer-Dewar.

#include "mfe.h"
#include "parametrizer_types.h"
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
        ("outfile,o", po::value<std::string>(), "Output structure file")
        ("paramdir,p", po::value<std::string>()->default_value("Turner99"), "Turner99 parameter directory")
        ("multiloop-penalty,a", po::value<std::string>(), "Multiloop penalty parameter")
        ("unpaired-penalty,b", po::value<std::string>(), "Unpaired base penalty parameter")
        ("branch-penalty,c", po::value<std::string>(), "Branching helix penalty parameter")
        ("dummy-scaling,d", po::value<std::string>(), "Dummy scaling parameter")
        ("help,h", "Display this help message")
        ;

    po::positional_options_description p;
    p.add("sequence", -1);
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    };

    // Process file-related options
    fs::path seq_file, output_file, param_dir;

    seq_file = fs::path(vm["sequence"].as<std::string>());
    param_dir = fs::path(vm["paramdir"].as<std::string>());

    if (vm.count("outfile")) {
        output_file = fs::path(vm["outputfile"].as<std::string>());
    } else {
        output_file = seq_file;
        output_file.replace_extension(".ct");
    }

    // Set up the parameter vector
    ParameterVector params = ParameterVector();

    if (vm.count("multiloop-penalty")) {
        params.multiloop_penalty = get_mpq_from_word(vm["multiloop-penalty"].as<std::string>());
    };

    if (vm.count("unpaired-penalty")) {
        params.unpaired_penalty = get_mpq_from_word(vm["unpaired-penalty"].as<std::string>());
    };

    if (vm.count("branch-penalty")) {
        params.branch_penalty = get_mpq_from_word(vm["branch-penalty"].as<std::string>());
    };

    if (vm.count("dummy-scaling")) {
        params.dummy_scaling = get_mpq_from_word(vm["dummy-scaling"].as<std::string>());
    };

    ScoreVector result = mfe(seq_file.native(), output_file.native(), param_dir.native(), params);
    std::cout << result << std::endl;
    return(0);
}
