// Copyright (c) 2014 Andrew Gainer-Dewar.

#include "nndb_constants.h"
#include "nntm.h"
#include "pmfe_types.h"
#include "rna_polytope.h"
#include "thread_pool.h"

#include <string>

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

#include <boost/mpi.hpp>

#define BOOST_LOG_DYN_LINK 1 // Fix an issue with dynamic library loading
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace mpi = boost::mpi;

int main(int argc, char * argv[]) {
    // Set up mpi
    mpi::environment env;
    mpi::communicator world;

    // Set up argument processing
    po::options_description desc("Options");
    desc.add_options()
        ("sequence", po::value<std::string>()->required(), "Sequence file")
        ("verbose,v", po::bool_switch()->default_value(false), "Write verbose debugging output")
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

    // Process logging-related options
    bool verbose = vm["verbose"].as<bool>();
    if (verbose) {
        boost::log::core::get()->set_filter(
            boost::log::trivial::severity >= boost::log::trivial::info);
    } else {
        boost::log::core::get()->set_filter
            (boost::log::trivial::severity >= boost::log::trivial::warning);
    }

    // Set up dangle model
    pmfe::dangle_mode dangles = pmfe::convert_to_dangle_mode(vm["dangle-model"].as<int>());

    // Set up sequence
    fs::path seq_file (vm["sequence"].as<std::string>());
    pmfe::RNASequence sequence(seq_file);

    if (world.rank() == 0) { // Master node
        // Build the polytope
        pmfe::RNAPolytope poly(sequence, dangles, thread_pool);

        poly.build(world.size() - 1);

        // Output the results
        poly.print_statistics();

        fs::path poly_file;
        if (vm.count("outfile")) {
            poly_file = fs::path(vm["outfile"].as<std::string>());
        } else {
            poly_file = seq_file;
            poly_file.replace_extension(".rnapoly");
        }

        poly.write_to_file(poly_file);
    } else { // Worker node
        while (true) {
            // Receive the question
            pmfe::Question question;
            world.recv(0, pmfe::QUESTION, question);

            // Quit if it's quittin' time
            if (question.kill_signal == true)
                break;

            // Otherwise, set up the computational environment
            pmfe::ParameterVector params = question.as_pv();
            pmfe::Turner99 constants(thread_pool, params);
            pmfe::NNTM energy_model(constants, dangles, thread_pool);

            // Compute the energy tables
            pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(sequence);

            // Find the MFE structure
            pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

            // Construct and return the answer
            pmfe::Answer answer (scored_structure);
            world.send(0, pmfe::ANSWER, answer);
        }
    }

    return 0;
}
