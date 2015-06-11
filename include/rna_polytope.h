// Copyright (c) 2015 Andrew Gainer-Dewar.

#ifndef RNA_POLYTOPE_H
#define RNA_POLYTOPE_H

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
#include <utility>

#include <boost/mpi.hpp>

#include "pmfe_types.h"
#include "rna_polytope.h"
#include "thread_pool.h"

#include <map>
#include "boost/filesystem/fstream.hpp"

#include <CGAL/Gmpq.h>

namespace fs = boost::filesystem;

namespace pmfe {
    enum MPI_JOBS {QUESTION, ANSWER};
    typedef CGAL::Gmpq Q; // We'll do all the geometry over Q

    class RNAPolytope: public CGAL::Convex_hull_d< CGAL::Cartesian_d<Q> > {
    public:
        // Typedef magic, part 1
        typedef CGAL::Cartesian_d<Q> K;
        typedef CGAL::Convex_hull_d<K> ConvexHull;
        typedef CGAL::Point_d< CGAL::Cartesian_d<Q> > FPoint;
        typedef CGAL::Vector_d<K> FVector;
        typedef typename ConvexHull::R R;
        typedef typename ConvexHull::Hull_vertex_iterator Hull_vertex_iterator;

        RNAPolytope(RNASequence sequence, dangle_mode dangles, SimpleThreadPool& thread_pool);

        FPoint vertex_oracle(FVector objective);
        void build(unsigned int n_parallel_tests);
        void write_to_file(const fs::path poly_file) const;

    protected:
        class compare_fp {
        public:
            bool operator() (const FPoint& left, const FPoint& right) const {
                // Custom comparitor for the results map in the polytope
                for (int i = 0; i < left.dimension(); ++i) {
                    if (left[i] < right[i]) {
                        return true;
                    } else if (left[i] > right[i]) {
                        return false;
                    }
                    // Otherwise, continue
                }

                // If we make it out of the loop, the points are equal
                return false;
            }
        };

        // Typedef magic, part 2
        typedef CGAL::Hyperplane_d<K> Hyperplane;
        typedef CGAL::Linear_algebraCd<Q> LinearAlgebra;
        typedef typename LinearAlgebra::Matrix LMatrix;
        typedef typename LinearAlgebra::Vector LVector;
        typedef typename ConvexHull::Facet_iterator Facet_iterator;
        typedef typename ConvexHull::Facet_handle Facet_handle;
        typedef typename FVector::Base_vector Base_vector;

        ScoreVector classical_scores;
        RNASequence sequence;
        dangle_mode dangles;
        std::map<FPoint, RNAStructureWithScore, compare_fp> structures;

        SimpleThreadPool& thread_pool;

        void hook_preinit();
        void hook_postinit();
        void hook_perloop(unsigned int confirmed, unsigned int facets);
        void hook_postloop();
        void hook_unconfirmed(Facet_handle facet) {};
        void hook_confirmed(Facet_handle facet) {};

        std::pair<unsigned int, unsigned int> facet_counts();
    };

    class Question {
    public:
        Question() {};

    Question(ParameterVector params):
        multiloop_penalty(params.multiloop_penalty.get_str()),
            unpaired_penalty(params.unpaired_penalty.get_str()),
            branch_penalty(params.branch_penalty.get_str()),
            dummy_scaling(params.dummy_scaling.get_str())
            {};

    Question(bool kill_signal):
        kill_signal(kill_signal)
        {};

        std::string multiloop_penalty;
        std::string unpaired_penalty;
        std::string branch_penalty;
        std::string dummy_scaling;
        bool kill_signal = false;

        ParameterVector as_pv() const {
            mpq_class a(multiloop_penalty), b(unpaired_penalty), c(branch_penalty), d(dummy_scaling);
            ParameterVector result(a, b, c, d);
            return result;
        }

    private:
        friend class boost::serialization::access;
        template <class Archive>
            void serialize(Archive & ar, const unsigned int ve_sion) {
            ar & multiloop_penalty;
            ar & unpaired_penalty;
            ar & branch_penalty;
            ar & dummy_scaling;
            ar & kill_signal;
        }
    };

    class Answer {
    public:
        Answer() {};

    Answer(RNAStructureWithScore scored_structure):
        multiloops(scored_structure.score.multiloops.get_si()),
            branches(scored_structure.score.branches.get_si()),
            unpaired(scored_structure.score.unpaired.get_si()),
            w(scored_structure.score.w.get_str()),
            energy(scored_structure.score.energy.get_str()),
            structure_txt(scored_structure.string()),
            seq_txt(scored_structure.seq.string())
            {};

        unsigned int multiloops;
        unsigned int branches;
        unsigned int unpaired;
        std::string w;
        std::string energy;
        std::string structure_txt;
        std::string seq_txt;

        RNAStructureWithScore as_struct() const {
            RNAStructure structure (seq_txt, structure_txt);
            ScoreVector score (
                multiloops,
                unpaired,
                branches,
                mpq_class(w),
                mpq_class(energy)
                );
            RNAStructureWithScore result (structure, score);
            return result;
        }

    private:
        friend class boost::serialization::access;
        template <class Archive>
            void serialize(Archive & ar, const unsigned int version) {
            ar & multiloops;
            ar & branches;
            ar & unpaired;
            ar & w;
            ar & energy;
            ar & structure_txt;
            ar & seq_txt;
        }
    };
}
#endif
