// Copyright (c) 2014 Andrew Gainer-Dewar.

#include "BBPolytope.h"
#include "mfe.h"
#include "parametrizer_types.h"
#include <iostream>
#include <string>
#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

#include <CGAL/Gmpq.h>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

typedef CGAL::Gmpq Q; // We'll do all the geometry over Q
typedef iB4e::BBPolytope<Q> BBP;

ParameterVector fv_to_pv(BBP::FVector v) {
    ParameterVector pv (mpq_class(v.cartesian(0).mpq()), mpq_class(v.cartesian(1).mpq()), mpq_class(v.cartesian(2).mpq()), mpq_class(v.cartesian(3).mpq()));
    return pv;
};

BBP::FPoint sv_to_fp(ScoreVector v) {
    mpq_t values [4];
    mpq_inits(values[0], values[1], values[2], values[3], NULL);

    mpq_set_z(values[0], v.multiloops.get_mpz_t());
    mpq_set_z(values[1], v.unpaired.get_mpz_t());
    mpq_set_z(values[2], v.branches.get_mpz_t());
    mpq_set(values[3], v.w.get_mpq_t());

    BBP::FPoint result(4, values, values+4);
    mpq_clears(values[0], values[1], values[2], values[3], NULL);
    return result;
};


class RNAPolytope: public BBP {
    fs::path seq_file;
    fs::path param_dir;
    fs::path struct_dir;
    int dangle_model;

public:
    RNAPolytope(fs::path seq_file, fs::path param_dir, fs::path struct_dir, int dangle_model, int dim = 4):
        BBPolytope(dim),
        seq_file(seq_file),
        param_dir(param_dir),
        struct_dir(struct_dir),
        dangle_model(dangle_model)
    {};

    FPoint vertex_oracle(FVector objective) {
        std::string structure_ext = ".ct";
        fs::path initial_output_file = struct_dir / seq_file.stem();
        initial_output_file.replace_extension(structure_ext);

        ParameterVector params = fv_to_pv(objective);
        ScoreVector scores = mfe(seq_file.string(), initial_output_file.string(), params, param_dir.string(), dangle_model);

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

        return sv_to_fp(scores);
    };
};


int main(int argc, char * argv[]) {
    // Set up argument processing
    po::options_description desc("Options");
    desc.add_options()
        ("sequence", po::value<std::string>()->required(), "Sequence file")
        ("paramdir,p", po::value<std::string>()->default_value("Turner99"), "Turner99 parameter directory")
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

    // Process file-related options
    fs::path seq_file, param_dir;

    // Setup dangle model
    int dangle_model = vm["dangle-model"].as<int>();

    seq_file = fs::path(vm["sequence"].as<std::string>());
    param_dir = fs::path(vm["paramdir"].as<std::string>());

    fs::path struct_dir = seq_file.parent_path() / seq_file.stem();
    fs::remove_all(struct_dir);
    fs::create_directory(struct_dir);

    RNAPolytope *poly = new RNAPolytope(seq_file, param_dir, struct_dir, dangle_model);
    poly->build();

    poly->print_statistics();

    return 0;
};
