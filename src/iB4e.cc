// Copyright (c) 2014 Andrew Gainer-Dewar.

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <set>
#include <string>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/lexical_cast.hpp>

#include "iB4e.h"
#include "mfe.h"
#include "parametrizer_types.h"

#include "NLTemplate.h"

namespace fs = boost::filesystem;

QPoint run_gtmfe_on_q4vector(QVector param_vector, fs::path seq_file, fs::path output_dir, fs::path param_dir, int dangle_model) {
    ParameterVector params;
    ScoreVector scores;

    std::string structure_ext = ".ct";
    fs::path seq_base_name = seq_file.stem();

    fs::path initial_output_file = output_dir / seq_base_name;
    initial_output_file.replace_extension(structure_ext);

    params = ParameterVector(param_vector);
    scores = mfe(seq_file.string(), initial_output_file.string(), param_dir.string(), params, dangle_model);

    std::string score_sep (", ");
    std::string w_score_string (boost::lexical_cast<std::string>(mpf_class(scores.w).get_d()));
    std::string score_string =
        ".[" +
        scores.multiloops.get_str(10) + score_sep +
        scores.unpaired.get_str(10) + score_sep +
        scores.branches.get_str(10) + score_sep +
        w_score_string +
        "]";

    fs::path updated_output_file = initial_output_file;
    updated_output_file.replace_extension(score_string + structure_ext);
    fs::rename(initial_output_file, updated_output_file);

    return scores.get_q4point();
}

int iB4e_main(std::string seq_file_path, std::string param_dir_path, int dangle_model) {
    // Work in four dimensions
    int dim = 4;

    // Initialize the convex hull object
    ConvexHull CH(dim);

    // Create directory path objects
    fs::path param_dir(param_dir_path);

    fs::path seq_file(seq_file_path);
    fs::path seq_base_name = seq_file.stem();

    fs::path struct_dir = seq_file.parent_path() / seq_base_name;
    fs::remove_all(struct_dir);
    fs::create_directory(struct_dir);

    // Fancy functional trick to create a vertex oracle function
    boost::function<QPoint (QVector param_vector)> vertex_finder =
        boost::bind(&run_gtmfe_on_q4vector, _1, seq_file, struct_dir, param_dir, dangle_model);

    // BEGIN LOGIC
    // INITIALIZATION

    // We keep a list of vectors which have been tested or are in the polytope
    std::vector<LVector> test_vectors;

    // To bootstrap the process, we compute the classical scores here
    // (we'll need them later anyway)
    QVector classical = ParameterVector().as_QVector();
    QPoint classical_mfe = vertex_finder(classical);
    CH.insert(classical_mfe);
    LVector classical_structure_point = LVector(classical_mfe.cartesian_begin(), classical_mfe.cartesian_end());

    // Huggins' initialization loop
    while (CH.current_dimension() < dim) {
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
            QVector c(dim, x.begin(), x.end());
            LVector new_test_vector;

            // If either of c or -c yields a point that increases the dimension of CH, add it
            if (!done) {
                QPoint mfe_point = vertex_finder(c);
                if (CH.number_of_vertices() == 0 || CH.is_dimension_jump(mfe_point)) {
                    CH.insert(mfe_point);
                    new_test_vector = LVector(mfe_point.cartesian_begin(), mfe_point.cartesian_end()) - classical_structure_point;
                    done = true;
                }
            }

            if (!done) {
                QPoint mfe_point = vertex_finder(-c);
                if (CH.is_dimension_jump(mfe_point)) {
                    CH.insert(mfe_point);
                    new_test_vector = LVector(mfe_point.cartesian_begin(), mfe_point.cartesian_end()) - classical_structure_point;
                    done = true;
                }
            }

            // Otherwise, add c as a test vector
            if (!done) {
                new_test_vector = x;
                done = true;
            }

            test_vectors.push_back(new_test_vector);
        } else { // If no orthogonal vector was found, the polytope is not full-dimensional and the algorithm will fail
            std::cerr << "Failed to bootstrap full-dimensional polytope!" << std::endl;
            exit(1);
        }
    }

    assert(CH.is_valid());
    assert(CH.current_dimension() == dim);

    // MAIN LOOP
    bool all_confirmed_so_far = false;
    int confirmed = 0;
    while (!all_confirmed_so_far) {
        all_confirmed_so_far = true;

        // Attempt to confirm every facet
        for (Facet_iterator f = CH.facets_begin(); f != CH.facets_end() && all_confirmed_so_far ; f++) {
            if (!f->is_confirmed()) { // If the facet is not already confirmed, test it
                Hyperplane hp = CH.hyperplane_supporting(f);
                QVector innernormal = -hp.orthogonal_vector(); // CGAL returns the outer normal
                //innernormal *= innernormal[0].denominator();
                QPoint mfe_point = vertex_finder(innernormal);

                switch (hp.oriented_side(mfe_point)) {
                case CGAL::ON_POSITIVE_SIDE:
                    CH.insert(mfe_point); // If so, add the new point and break the inner loop
                    all_confirmed_so_far = false;
                    break;

                case CGAL::ON_ORIENTED_BOUNDARY:
                    f->confirm(); // Otherwise, mark the facet as confirmed and continue
                    confirmed++;
                    break;

                case CGAL::ON_NEGATIVE_SIDE:
                    std::cerr << "Failed vector test!" << std::endl;
                    std::cerr << "Hyperplane: " << hp << std::endl;
                    std::cerr << "Inner normal: " << innernormal << std::endl;
                    std::cerr << "MFE point: " << mfe_point << std::endl;
                    std::cerr << "MFE energy: " << std::inner_product(mfe_point.cartesian_begin(), mfe_point.cartesian_end(), innernormal.cartesian_begin(), CGAL::Gmpq(0)) << std::endl;
                    QPoint known_v =  CH.vertex_of_facet(f, 0)->point();
                    std::cerr << "Example known vertex: " << known_v << std::endl;
                    std::cerr << "Known vertex energy: " << std::inner_product(known_v.cartesian_begin(), known_v.cartesian_end(), innernormal.cartesian_begin(), CGAL::Gmpq(0)) << std::endl << std::endl;
                    f->confirm();
                    break;
                }
            }
        }
    }

    // END LOGIC
    // BEGIN OUTPUT

    // Output statistics to stdout
    assert(CH.is_valid());
    CH.print_statistics();
    std::cout << CH.all_facets().size() << " facets, " << confirmed << " confirmed." << std::endl;

    // Build string encoding discovered points in polytope
    std::ostringstream point_stream;
    point_stream << "[";

    for (Vertex_iterator v = CH.vertices_begin(); v != CH.vertices_end(); v++) {
        QPoint p = v->point();
        point_stream << "[" << p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << "], ";
    }

    point_stream << "]";

    // Build string encoding classical scores
    std::ostringstream classical_score_stream;
    classical_score_stream << "[" <<
        classical_structure_point[0] << ", " <<
        classical_structure_point[1] << ", " <<
        classical_structure_point[2] << ", " <<
        classical_structure_point[3] << "]";

    // Build polytope file from template
    NL::Template::LoaderFile loader;
    NL::Template::Template t(loader);

    t.load("output.template");
    t.set("POINTS", point_stream.str());
    t.set("CLASSICAL_SCORES", classical_score_stream.str());
    t.set("POLYDUMP", seq_base_name.string());

    NL::Template::OutputString result_buffer;
    t.render(result_buffer);

    fs::path polytope_file_path = seq_file;
    polytope_file_path.replace_extension(".polytope.sage");

    fs::ofstream polytope_file(polytope_file_path);
    polytope_file << result_buffer.buf.str();

    // END OUTPUT

    return 0;
}
