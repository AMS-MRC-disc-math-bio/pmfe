// Copyright (c) 2014 Andrew Gainer-Dewar.

#ifndef BBPOLY_H
#define BBPOLY_H

#include <CGAL/Cartesian_d.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/Regular_complex_d.h>

#include <CGAL/Kernel_d/Point_d.h>
#include <CGAL/Kernel_d/Vector_d.h>
#include <CGAL/Kernel_d/Hyperplane_d.h>
#include <CGAL/Linear_algebraCd.h>

#include <CGAL/Origin.h>
#include <CGAL/Gmpq.h>

#include <string>

namespace iB4e
{

    template <typename F>
        class BBPolytope: public CGAL::Convex_hull_d< CGAL::Cartesian_d<F> > {

    public:
        // Typedef magic, part 1
        typedef CGAL::Cartesian_d<F> K;
        typedef CGAL::Convex_hull_d<K> ConvexHull;
        typedef CGAL::Point_d<K> FPoint;
        typedef CGAL::Vector_d<K> FVector;
        typedef typename ConvexHull::R R;
        typedef typename ConvexHull::Hull_vertex_iterator Hull_vertex_iterator;

    BBPolytope(int dim):
        ConvexHull(dim, R())
        {};

        void build(); // Implementation below for readability
        virtual FPoint vertex_oracle (FVector objective) = 0;

    protected:
        // Typedef magic, part 2
        typedef CGAL::Hyperplane_d<K> Hyperplane;
        typedef CGAL::Linear_algebraCd<F> LinearAlgebra;
        typedef typename LinearAlgebra::Matrix LMatrix;
        typedef typename LinearAlgebra::Vector LVector;
        typedef typename ConvexHull::Facet_iterator Facet_iterator;
        typedef typename FVector::Base_vector Base_vector;
    };

    template <typename F>
        void BBPolytope<F>::build() {
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

        // Huggins' initialization loop
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
                if (!done) {
                    FPoint result = this->vertex_oracle(c);
                    if (this->number_of_vertices() == 0 || this->is_dimension_jump(result)) {
                        this->insert(result);
                        new_test_vector = LVector(result.cartesian_begin(), result.cartesian_end()) - first_result_lv;
                        done = true;
                    }
                }

                if (!done) {
                    FPoint result = vertex_oracle(-c);
                    if (this->is_dimension_jump(result)) {
                        this->insert(result);
                        new_test_vector = LVector(result.cartesian_begin(), result.cartesian_end()) - first_result_lv;
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
                // TODO: Throw an exception instead
                exit(-1);
            }
        }

        assert(this->is_valid());
        assert(this->current_dimension() == dim);

        // MAIN LOOP
        bool all_confirmed_so_far = false;
        int confirmed = 0;
        while (!all_confirmed_so_far) {
            all_confirmed_so_far = true;

            // Attempt to confirm every facet
            for (Facet_iterator f = this->facets_begin(); f != this->facets_end() && all_confirmed_so_far ; f++) {
                if (!f->is_confirmed()) { // If the facet is not already confirmed, test it
                    Hyperplane hp = this->hyperplane_supporting(f);
                    FVector innernormal = -hp.orthogonal_vector(); // CGAL returns the outer normal
                    //innernormal *= innernormal[0].denominator();
                    FPoint result = vertex_oracle(innernormal);

                    switch (hp.oriented_side(result)) {
                    case CGAL::ON_POSITIVE_SIDE:
                        this->insert(result); // If so, add the new point and break the inner loop
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
                        std::cerr << "Found point: " << result << std::endl;
                        std::cerr << "Found point score: " << std::inner_product(result.cartesian_begin(), result.cartesian_end(), innernormal.cartesian_begin(), CGAL::Gmpq(0)) << std::endl;
                        FPoint known_v =  this->vertex_of_facet(f, 0)->point();
                        std::cerr << "Example known vertex: " << known_v << std::endl;
                        std::cerr << "Known vertex score: " << std::inner_product(known_v.cartesian_begin(), known_v.cartesian_end(), innernormal.cartesian_begin(), CGAL::Gmpq(0)) << std::endl << std::endl;
                        f->confirm();
                        break;
                    }
                }
            }
        }
        // END LOGIC
    };
}
#endif
