// Copyright (c) 2015 Andrew Gainer-Dewar.

#include "rna_polytope.h"
#include "pmfe_types.h"
#include "nntm.h"
#include "nndb_constants.h"
#include "BBPolytope.h"
#include "mfe.h"
#include "thread_pool.h"

#include <map>

#include <CGAL/Gmpq.h>

#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"

#define BOOST_LOG_DYN_LINK 1 // Fix an issue with dynamic library loading
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

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
        std::vector<Rational> values = {structure.score.multiloops, structure.score.unpaired, structure.score.branches, structure.score.w};
        BBP::FPoint result(4, values.begin(), values.end());
        return result;
    };

    RNAPolytope::RNAPolytope(RNASequence sequence, pmfe::dangle_mode dangles, SimpleThreadPool& thread_pool):
        BBPolytope(4),
        sequence(sequence),
        dangles(dangles),
        thread_pool(thread_pool)
        {};

    BBP::FPoint RNAPolytope::vertex_oracle(FVector objective) {
        // Set up the computational environment
        ParameterVector params = fv_to_pv(objective);
        Turner99 constants(thread_pool, params);
        NNTM energy_model(constants, dangles, thread_pool);

        // Compute the energy tables
        RNASequenceWithTables seq_annotated = energy_model.energy_tables(sequence);

        // Find the MFE structure
        RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);
        BBP::FPoint result = scored_structure_to_fp(scored_structure);

        // TODO: Handle storing stuctures in class after conversion to dD_triangulation
        structures.insert(std::make_pair(result, scored_structure));
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

        outfile << "#\t" << sequence << "\tm\tu\th\tw\te" << std::endl;

        // Initializing variables outside the loop is bad manners, but
        // we need it in order to have both indices
        int i;
        BBP::Hull_vertex_const_iterator v;
        for (i = 1, v = hull_vertices_begin(); v != hull_vertices_end(); ++i, ++v) {
            outfile << i << "\t" << structures.at(associated_point(v)) << std::endl;
        }
    };

    void RNAPolytope::hook_preinit() {
        BOOST_LOG_TRIVIAL(info) << "Initializing polytope.";
    };

    void RNAPolytope::hook_postinit() {
        BOOST_LOG_TRIVIAL(info) << "Initialization complete. Beginning loop.";
    };

    void RNAPolytope::hook_perloop(size_t confirmed) {
        BOOST_LOG_TRIVIAL(info) << confirmed << " / " << number_of_simplices() << " known facets confirmed.";
    };

    void RNAPolytope::hook_postloop() {
        BOOST_LOG_TRIVIAL(info) << "Polytope complete.";
    };
};
