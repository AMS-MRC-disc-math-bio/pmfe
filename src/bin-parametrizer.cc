// Copyright (c) 2014 Andrew Gainer-Dewar.

#include "BBPolytope.h"
#include "mfe.h"
#include "nndb_constants.h"
#include "nntm.h"
#include "pmfe_types.h"

#include <iostream>
#include <string>
#include <fstream>

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

#include <CGAL/Gmpq.h>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

typedef CGAL::Gmpq Q; // We'll do all the geometry over Q
typedef iB4e::BBPolytope<Q, pmfe::RNAStructureWithScore> BBP;

pmfe::ParameterVector fv_to_pv(BBP::FVector v) {
    pmfe::ParameterVector pv (
        mpq_class(v.cartesian(0).mpq()),
        mpq_class(v.cartesian(1).mpq()),
        mpq_class(v.cartesian(2).mpq()),
        mpq_class(v.cartesian(3).mpq())
        );
    return pv;
};

BBP::FPoint scored_structure_to_fp(pmfe::RNAStructureWithScore structure) {
    mpq_t values [4];
    mpq_inits(values[0], values[1], values[2], values[3], NULL);

    mpq_set_z(values[0], structure.score.multiloops.get_mpz_t());
    mpq_set_z(values[1], structure.score.unpaired.get_mpz_t());
    mpq_set_z(values[2], structure.score.branches.get_mpz_t());
    mpq_set(values[3], structure.score.w.get_mpq_t());

    BBP::FPoint result(4, values, values+4);
    mpq_clears(values[0], values[1], values[2], values[3], NULL);
    result.payload = structure;
    return result;
};

class RNAPolytope: public BBP {
public:
    pmfe::ScoreVector classical_scores;
    pmfe::RNASequence sequence;
    pmfe::dangle_mode dangles;

    RNAPolytope(pmfe::RNASequence sequence, pmfe::dangle_mode dangles):
        BBPolytope(4),
        sequence(sequence),
        dangles(dangles)
        {};

    FPoint vertex_oracle(FVector objective) {
        // Set up the computational environment
        pmfe::ParameterVector params = fv_to_pv(objective);
        pmfe::Turner99 constants(params);
        pmfe::NNTM energy_model(constants, dangles);

        // Compute the energy tables
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(sequence);

        // Find the MFE structure
        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        // Compute the classical score of the structure
        pmfe::Turner99 classical_constants;
        pmfe::NNTM classical_model(classical_constants, dangles);
        mpq_class classical_energy = classical_model.energy(scored_structure);

        // Build the score vector
        scored_structure.score.w = classical_energy - (scored_structure.score.multiloops * classical_constants.multConst[0] + scored_structure.score.unpaired * classical_constants.multConst[1] + scored_structure.score.branches * classical_constants.multConst[2]);
        scored_structure.score.canonicalize();

        // Check that the w calculation produced a consistent result
        mpq_class formula_energy = scored_structure.score.multiloops * params.multiloop_penalty + scored_structure.score.unpaired * params.unpaired_penalty + scored_structure.score.branches * params.branch_penalty + scored_structure.score.w * params.dummy_scaling;
        formula_energy.canonicalize();

        // And alert the user if not
        if (scored_structure.score.energy != formula_energy) {
            std::cerr << "Energy calculation is inconsistent!" << std::endl;
            std::cerr << params << std::endl;
            std::cerr << scored_structure.score << std::endl;
            std::cerr << "Formula energy: " << formula_energy.get_str(10) << std::endl;
            std::cerr << "Classical energy: " << classical_energy.get_str(10) << std::endl << std::endl;
        };

        return scored_structure_to_fp(scored_structure);
    };
};

void write_poly_file(RNAPolytope* poly, const fs::path poly_file) {
    // TODO: Write this!
    std::cout << "Request to write polytope file acknowledged." << std::endl;
};


int main(int argc, char * argv[]) {
    // Set up argument processing
    po::options_description desc("Options");
    desc.add_options()
        ("sequence", po::value<std::string>()->required(), "Sequence file")
        ("outfile,o", po::value<std::string>(), "Output file")
        ("dangle-model,m", po::value<int>()->default_value(1), "Dangle model")
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

    // Set up dangle model
    pmfe::dangle_mode dangles = pmfe::convert_to_dangle_mode(vm["dangle-model"].as<int>());

    // Set up sequence
    fs::path seq_file (vm["sequence"].as<std::string>());
    pmfe::RNASequence sequence(seq_file);

    RNAPolytope* poly = new RNAPolytope(sequence, dangles);
    poly->build();

    poly->print_statistics();

    fs::path poly_file;
    if (vm.count("outfile")) {
        poly_file = fs::path(vm["outfile"].as<std::string>());
    } else {
        poly_file = seq_file;
        poly_file.replace_extension(".poly");
    }

    write_poly_file(poly, poly_file);

    return 0;
};
