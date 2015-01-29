// Copyright (c) 2014 Andrew Gainer-Dewar.

#include "BBPolytope.h"
#include "mfe.h"
#include "pmfe_types.h"
#include <iostream>
#include <string>
#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"
#include "NLTemplate.h"
#include <fstream>

#include <CGAL/Gmpq.h>

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace nlt = NL::Template;

typedef CGAL::Gmpq Q; // We'll do all the geometry over Q
typedef iB4e::BBPolytope<Q> BBP;

pmfe::ParameterVector fv_to_pv(BBP::FVector v) {
    pmfe::ParameterVector pv (
                              mpq_class(v.cartesian(0).mpq()),
                              mpq_class(v.cartesian(1).mpq()),
                              mpq_class(v.cartesian(2).mpq()),
                              mpq_class(v.cartesian(3).mpq())
                              );
    return pv;
};

BBP::FVector pv_to_fv(pmfe::ParameterVector v) {
    mpq_t values [4];
    mpq_inits(values[0], values[1], values[2], values[3], NULL);

    mpq_set(values[0], v.multiloop_penalty.get_mpq_t());
    mpq_set(values[1], v.unpaired_penalty.get_mpq_t());
    mpq_set(values[2], v.branch_penalty.get_mpq_t());
    mpq_set(values[3], v.dummy_scaling.get_mpq_t());

    BBP::FVector result(4, values, values+4);
    mpq_clears(values[0], values[1], values[2], values[3], NULL);
    return result;
};

pmfe::ScoreVector fp_to_sv(BBP::FPoint v) {
    pmfe::ScoreVector sv (
                          mpz_class(mpq_class(v.cartesian(0).mpq())),
                          mpz_class(mpq_class(v.cartesian(1).mpq())),
                          mpz_class(mpq_class(v.cartesian(2).mpq())),
                          mpq_class(v.cartesian(3).mpq())
                          );
    return sv;
};

BBP::FPoint sv_to_fp(pmfe::ScoreVector v) {
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
    pmfe::ScoreVector classical_scores;

    RNAPolytope(fs::path seq_file, fs::path param_dir, fs::path struct_dir, int dangle_model, int dim = 4):
        BBPolytope(dim),
        seq_file(seq_file),
        param_dir(param_dir),
        struct_dir(struct_dir),
        dangle_model(dangle_model)
    {
        init();
    };

    FPoint vertex_oracle(FVector objective) {
        std::string structure_ext = ".ct";
        fs::path initial_output_file = struct_dir / seq_file.stem();
        initial_output_file.replace_extension(structure_ext);

        pmfe::ParameterVector params = fv_to_pv(objective);
        pmfe::ScoreVector scores = pmfe::mfe(seq_file.string(), initial_output_file.string(), params, param_dir.string(), dangle_model);

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

private:
    void init() {
        pmfe::ParameterVector classical_params;
        fs::path classical_output_file = struct_dir / seq_file.stem();
        classical_output_file.replace_extension(".classical.ct");
        classical_scores = pmfe::mfe(seq_file.string(), classical_output_file.string(), classical_params, param_dir.string(), dangle_model);
    };
};

void create_poly_file(RNAPolytope* poly, fs::path template_file, fs::path poly_file) {
    // Load the output template
    nlt::LoaderFile loader;
    nlt::Template t(loader);
    t.load(template_file.c_str());

    // Set relevant variables
    fs::path poly_base = poly_file.stem();
    t.set("SEQNAME", poly_base.string());
    t.set("POLYDUMP", poly_base.string());
    t.set("POINTCOUNT", boost::lexical_cast<std::string>(poly->number_of_vertices()));
    t.set("FACETCOUNT", boost::lexical_cast<std::string>(poly->number_of_simplices()));
    t.set("CLASSICAL_SCORES", poly->classical_scores.print_as_list());

    // Build the Python list of points
    std::string pointslist = "[";
    for (BBP::Hull_vertex_iterator v = poly->hull_vertices_begin(); v != poly->hull_vertices_end(); ++v) {
        pmfe::ScoreVector sv = fp_to_sv(v->point());
        pointslist += sv.print_as_list();
        pointslist += ", ";
    }
    pointslist += "]";
    t.set("POINTS", pointslist);

    // Render the results to an NLTemplate OutputString
    nlt::OutputString result;
    t.render(result);

    // Write the string to the output file
    std::ofstream outfile (poly_file.c_str());
    outfile << result.buf.str();
    outfile.close();

    // Delete the sobj file if it exists (so we don't accidentally use an outdated one)
    fs::path poly_sobj = poly_base;
    poly_sobj += ".sobj";
    fs::remove(poly_sobj);
};


int main(int argc, char * argv[]) {
    // Set up argument processing
    po::options_description desc("Options");
    desc.add_options()
        ("sequence", po::value<std::string>()->required(), "Sequence file")
        ("paramdir,p", po::value<std::string>()->default_value("/usr/local/share/pmfe/Turner99/pmfe"), "Turner99 parameter directory")
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
    fs::path template_file = "output.template";

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

    fs::path poly_file = seq_file.replace_extension(".sage");
    create_poly_file(poly, template_file, poly_file);

    return 0;
};
