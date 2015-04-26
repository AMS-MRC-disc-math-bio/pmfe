// Copyright (c) 2015 Andrew Gainer-Dewar.

#include "rna_polytope.h"
#include "pmfe_types.h"
#include "nntm.h"
#include "nndb_constants.h"
#include "BBPolytope.h"
#include "mfe.h"

#include <map>

#include <CGAL/Gmpq.h>

#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"

namespace fs = boost::filesystem;

namespace pmfe {
    ParameterVector fv_to_pv(BBP::FVector v) {
        ParameterVector pv (
            mpq_class(v.cartesian(0).mpq()),
            mpq_class(v.cartesian(1).mpq()),
            mpq_class(v.cartesian(2).mpq()),
            mpq_class(v.cartesian(3).mpq())
            );
        return pv;
    };

    BBP::FPoint scored_structure_to_fp(RNAStructureWithScore structure) {
        mpq_t values [4];
        mpq_inits(values[0], values[1], values[2], values[3], NULL);

        mpq_set_z(values[0], structure.score.multiloops.get_mpz_t());
        mpq_set_z(values[1], structure.score.unpaired.get_mpz_t());
        mpq_set_z(values[2], structure.score.branches.get_mpz_t());
        mpq_set(values[3], structure.score.w.get_mpq_t());

        BBP::FPoint result(4, values, values+4);
        mpq_clears(values[0], values[1], values[2], values[3], NULL);
        return result;
    };

    RNAPolytope::RNAPolytope(RNASequence sequence, pmfe::dangle_mode dangles):
        BBPolytope(4),
        sequence(sequence),
        dangles(dangles)
        {};

    BBP::FPoint RNAPolytope::vertex_oracle(FVector objective) {
        // Set up the computational environment
        ParameterVector params = fv_to_pv(objective);
        Turner99 constants(params);
        NNTM energy_model(constants, dangles);

        // Compute the energy tables
        RNASequenceWithTables seq_annotated = energy_model.energy_tables(sequence);

        // Find the MFE structure
        RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        // Compute the classical score of the structure
        Turner99 classical_constants;
        NNTM classical_model(classical_constants, dangles);
        ScoreVector classical_score = classical_model.score(scored_structure);
        mpq_class classical_energy = classical_score.energy;

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
            std::cerr << scored_structure << std::endl;
            std::cerr << "Formula energy: " << formula_energy.get_str(10) << std::endl;
            std::cerr << "Score energy: " << scored_structure.score.energy.get_str(10) << std::endl;
        };

        BBP::FPoint result = scored_structure_to_fp(scored_structure);
        // TODO: Handle storing stuctures in class after conversion to dD_triangulation
        structures[result] = scored_structure;
        return result;
    };

    void RNAPolytope::write_to_file(const fs::path poly_file) const {
        fs::ofstream outfile(poly_file);

        if(!outfile.is_open()) {
            std::cerr << "Couldn't open polytope file " << poly_file << "." << std::endl;
            exit(EXIT_FAILURE);
        }

        outfile << "# Points: " << number_of_vertices() << std::endl;
        outfile << "# Facets: " << number_of_simplices() << std::endl << std::endl;

        outfile << "#\t" << sequence << "m\tu\th\tw\te" << std::endl;

        // Initializing variables outside the loop is bad manners, but
        // we need it in order to have both indices
        int i;
        BBP::Hull_vertex_const_iterator v;
        for (i = 1, v = hull_vertices_begin(); v != hull_vertices_end(); ++i, ++v) {
            outfile << i << "\t" << structures.at(associated_point(v)) << std::endl;
        }
    };
};
