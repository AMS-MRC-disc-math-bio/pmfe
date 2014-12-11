// Copyright (c) 2014 Andrew Gainer-Dewar.

#include "iB4e.h"
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
        ("paramdir,p", po::value<std::string>()->default_value("Turner99"), "Turner99 parameter directory")
        ("dangle-model,m", po::value<int>(), "Dangle model")
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

    // Process file-related options
    fs::path seq_file, param_dir;

    // Setup dangle model
    int dangle_model = 1;
    if (vm.count("dangle-model")) {
        dangle_model = vm["dangle-model"].as<int>();;
    }

    seq_file = fs::path(vm["sequence"].as<std::string>());
    param_dir = fs::path(vm["paramdir"].as<std::string>());

    return iB4e_main(seq_file.native(), param_dir.native(), dangle_model);
}
