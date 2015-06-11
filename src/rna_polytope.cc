// Copyright (c) 2015 Andrew Gainer-Dewar.

#include "rna_polytope.h"
#include "pmfe_types.h"
#include "nntm.h"
#include "nndb_constants.h"
#include "mfe.h"
#include "thread_pool.h"

#include <map>
#include <cassert>
#include <utility>

#include <CGAL/Gmpq.h>

#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"

#include <boost/mpi.hpp>

#define BOOST_LOG_DYN_LINK 1 // Fix an issue with dynamic library loading
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

namespace fs = boost::filesystem;
namespace mpi = boost::mpi;

namespace pmfe {
    mpi::environment env;
    mpi::communicator world;

    ParameterVector fv_to_pv(RNAPolytope::FVector v) {
        ParameterVector pv (
            mpq_class(v.cartesian(0).mpq()),
            mpq_class(v.cartesian(1).mpq()),
            mpq_class(v.cartesian(2).mpq()),
            mpq_class(v.cartesian(3).mpq())
            );
        return pv;
    };

    RNAPolytope::FPoint scored_structure_to_fp(RNAStructureWithScore structure) {
        std::vector<Rational> values = {structure.score.multiloops, structure.score.unpaired, structure.score.branches, structure.score.w};
        RNAPolytope::FPoint result(4, values.begin(), values.end());
        return result;
    };

    RNAPolytope::RNAPolytope(RNASequence sequence, pmfe::dangle_mode dangles, SimpleThreadPool& thread_pool):
        ConvexHull(4, R()),
        sequence(sequence),
        dangles(dangles),
        thread_pool(thread_pool)
    {};

    RNAPolytope::FPoint RNAPolytope::vertex_oracle(RNAPolytope::FVector objective) {
        // Set up the computational environment
        ParameterVector params = fv_to_pv(objective);
        Turner99 constants(thread_pool, params);
        NNTM energy_model(constants, dangles, thread_pool);

        // Compute the energy tables
        RNASequenceWithTables seq_annotated = energy_model.energy_tables(sequence);

        // Find the MFE structure
        RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);
        RNAPolytope::FPoint result = scored_structure_to_fp(scored_structure);

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
        Hull_vertex_const_iterator v;
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

    void RNAPolytope::hook_perloop(unsigned int confirmed, unsigned int facets) {
        BOOST_LOG_TRIVIAL(info) << "Facets (confirmed / known): " << confirmed << " / " << facets << ".";
    };

    void RNAPolytope::hook_postloop() {
        std::pair<unsigned int, unsigned int> facet_counts = this->facet_counts();
        BOOST_LOG_TRIVIAL(info) << "Polytope complete: " << facet_counts.first << " / " << facet_counts.second;

        // Kill the workers
        pmfe::Question kill(true);

        for (unsigned int i = 1; i < world.size(); ++i) {
            world.send(i, QUESTION, kill);
        };
    };

    void RNAPolytope::build(unsigned int n_parallel_tests) {
        // Input specification
        assert (n_parallel_tests != 0);
        assert (n_parallel_tests <= world.size());

        // For now, only allow this to run if the polytope is empty
        //assert(this->current_dimension() == 0);
        int dim = this->dimension();

        // BEGIN LOGIC
        // INITIALIZATION
        // We'll keep a list of vectors which have been tested or are in the polytope
        std::vector<LVector> test_vectors;

        // Bootstrap the process by finding one point manually
        FVector first_test (dim, Base_vector(), 0);
        FPoint first_result = this->vertex_oracle(first_test);
        LVector first_result_lv = LVector(first_result.cartesian_begin(), first_result.cartesian_end());
        this->insert(first_result);

        // TODO: Why is this trying certain vectors many times?

        // Huggins' initialization loop
        hook_preinit();

        while (this->current_dimension() < dim) {
            LMatrix A;
            // Construct some helper objects
            if (test_vectors.size() > 0) {
                A = LMatrix(test_vectors);
            } else {
                A = LMatrix(dim, dim);
            }

            A = LinearAlgebra::transpose(A); // Transpose due to design choices in CGAL
            LMatrix orthogonal_vectors;
            LVector x = LVector(dim);
            bool done = false;

            // Find a vector orthogonal to all the test vectors
            if (LinearAlgebra::homogeneous_linear_solver(A, orthogonal_vectors)) {
                x = orthogonal_vectors.column(0);
                FVector c(dim, x.begin(), x.end());
                LVector new_test_vector;

                // If either of c or -c yields a point that increases the dimension of the polytope, add it
                if (not done) {
                    FPoint result = this->vertex_oracle(c);
                    if (this->number_of_vertices() == 0 or this->is_dimension_jump(result)) {
                        this->insert(result);
                        new_test_vector = LVector(result.cartesian_begin(), result.cartesian_end()) - first_result_lv;
                        done = true;
                    }
                }

                if (not done) {
                    FPoint result = vertex_oracle(-c);
                    if (this->is_dimension_jump(result)) {
                        this->insert(result);
                        new_test_vector = LVector(result.cartesian_begin(), result.cartesian_end()) - first_result_lv;
                        done = true;
                    }
                }

                // Otherwise, add c as a test vector
                if (not done) {
                    new_test_vector = x;
                    done = true;
                }

                test_vectors.push_back(new_test_vector);
            } else { // If no orthogonal vector was found, the polytope is not full-dimensional and the algorithm will fail
                std::cerr << "Failed to bootstrap full-dimensional polytope!" << std::endl;
                // TODO: Throw an exception instead
                exit(-1);
            }
        }

        assert(this->is_valid());
        assert(this->current_dimension() == dim);

        hook_postinit();

        // BEGIN MAIN LOOP
        while (true) {
            std::pair<unsigned int, unsigned int> facet_counts = this->facet_counts();
            unsigned int confirmed = facet_counts.first;
            unsigned int total_facets = facet_counts.second;
            if (confirmed == total_facets) {
                break;
            }

            hook_perloop(confirmed, total_facets);

            std::vector<Question> questions;
            std::vector<Facet_handle> facets;

            // Process up to n_parallel_tests facets
            for (Facet_iterator f = this->facets_begin(); f != this->facets_end() and questions.size() < n_parallel_tests; ++f) {
                if (not f->is_confirmed()) { // If the facet is not already confirmed, test it
                    facets.push_back(f);
                    Hyperplane hp = this->hyperplane_supporting(f);
                    FVector innernormal = -hp.orthogonal_vector(); // CGAL returns the outer normal

                    ParameterVector params = fv_to_pv(innernormal);
                    Question question(params);
                    questions.push_back(question);
                }
            }

            assert (questions.size() == facets.size());

            std::vector<mpi::request> sends;
            std::vector<mpi::request> recvs;
            std::vector<Answer> answers (questions.size());

            for (unsigned int i = 0; i < questions.size(); ++i) {
                sends.push_back(world.isend(i+1, QUESTION, questions[i]));
                recvs.push_back(world.irecv(i+1, ANSWER, answers[i]));
            };

            mpi::wait_all(sends.begin(), sends.end()); // Actually send the processing requests
            mpi::wait_all(recvs.begin(), recvs.end()); // Wait for the results

            assert (questions.size() == answers.size());

            // Once the results are obtained, we process them all
            for (unsigned int i = 0; i < questions.size(); ++i) {
                Facet_handle f = facets[i];
                Hyperplane hp = this->hyperplane_supporting(f);
                FVector innernormal = -hp.orthogonal_vector();

                RNAStructureWithScore scored_structure = answers[i].as_struct();
                FPoint result = scored_structure_to_fp(scored_structure);

                // Determine whether this result confirms or extends its facet
                switch (hp.oriented_side(result)) {
                case CGAL::ON_POSITIVE_SIDE: // The new point is outside the polytope
                    hook_unconfirmed(f);
                    this->insert(result); // so we add it as a vertex.
                    structures.insert(std::make_pair(result, scored_structure));
                    break;

                case CGAL::ON_ORIENTED_BOUNDARY: // The new point is on the boundary
                    hook_confirmed(f);
                    f->confirm(); // so the facet is confirmed.
                    structures.insert(std::make_pair(result, scored_structure));
                    break;

                case CGAL::ON_NEGATIVE_SIDE: // The new point is inside the polytope
                    // This is a logic error, so we make some noise
                    std::cerr << "Failed vector test!" << std::endl;
                    std::cerr << "Hyperplane: " << hp << std::endl;
                    std::cerr << "Inner normal: " << innernormal << std::endl;
                    std::cerr << "Found point: " << result << std::endl;
                    std::cerr << "Found point score: " << std::inner_product(result.cartesian_begin(), result.cartesian_end(), innernormal.cartesian_begin(), CGAL::Gmpq(0)) << std::endl;
                    FPoint known_v =  this->vertex_of_facet(f, 0)->point();
                    std::cerr << "Example known vertex: " << known_v << std::endl;
                    std::cerr << "Known vertex score: " << std::inner_product(known_v.cartesian_begin(), known_v.cartesian_end(), innernormal.cartesian_begin(), CGAL::Gmpq(0)) << std::endl << std::endl;
                    break;
                }
            };
        }
        // END MAIN LOOP

        hook_postloop();
        // END LOGIC
    };

    std::pair<unsigned int, unsigned int> RNAPolytope::facet_counts() {
        // First: number of confirmed facets
        // Second: total number of facets
        unsigned int n_confirmed = 0;
        unsigned int n_facets = 0;
        for (Facet_iterator f = this->facets_begin(); f != this->facets_end(); ++f) {
            ++n_facets;
            if (f->is_confirmed()) {
                ++n_confirmed;
            }
        }

        return std::make_pair(n_confirmed, n_facets);
    }
};
