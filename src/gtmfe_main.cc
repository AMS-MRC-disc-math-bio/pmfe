// Copyright (c) 2014 Andrew Gainer-Dewar.

#include "mfe.h"
#include "parametrizer_types.h"
#include <iostream>
#include <string>
#include "boost/program_options.hpp"
#include <gmpxx.h>

int main(int argc, char * argv[]) {
    namespace po = boost::program_options;
    po::options_description desc("Options");
    desc.add_options()
        ("sequence,s", po::value<std::string>(), "Sequence file")
        ("outfile,o", po::value<std::string>(), "Output structure file")
        ("paramdir,d", po::value<std::string>()->default_value("Turner99"), "Turner99 parameter directory")
        ("multiloop-penalty,a", po::value<std::string>(), "Multiloop penalty parameter")
        ("unpaired-penalty,b", po::value<std::string>(), "Unpaired base penalty parameter")
        ("branch-penalty,c", po::value<std::string>(), "Branching helix penalty parameter")
        ("dummy-scaling,d", po::value<std::string>(), "Dummy scaling parameter")
        ("help,h", "Display this help message")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    };

    ParameterVector params = ParameterVector();

    if (vm.count("multiloop-penalty")) {
        params.multiloop_penalty = mpq_class(vm["multiloop-penalty"].as<std::string>());
    };

    std::cout << mfe("combinatorial_seq.fasta", "combinatorial_seq.ct", "Turner99") << std::endl;
    return(0);
}
