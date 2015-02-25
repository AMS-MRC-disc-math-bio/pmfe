#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmpxx.h>

#include "loader.h"
#include "algorithms.h"
#include "mfe.h"
#include "subopt_traceback.h"
#include "global.h"
#include "utils.h"
#include "pmfe_types.h"

#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

namespace pmfe {
    using namespace std;

    static void write_header_subopt_file(const std::string subopt_file, const std::string& seq, const mpq_class energy)
    {
        ofstream outfile;
        outfile.open(subopt_file.c_str());
        char buff[4096];

        sprintf(buff,"count\t%s\t%f", seq.c_str(), energy.get_d());
        outfile << buff << std::endl;

        outfile.close();
    }

    static std::string subopt_filename(const std::string seq_file) {
        fs::path subopt_file = seq_file;

        subopt_file.replace_extension("ss.txt");

        return subopt_file.native();
    };

    void generate_subopt(const std::string seq_file, const std::string output_file, const mpq_class delta, const ParameterVector params, const std::string param_dir, const int max_structure_count, const int dangle_model) {
        std::string seq = "";
        if (read_sequence_file(seq_file.c_str(), seq) == FAILURE) {
            printf("Failed to open sequence file: %s.\n\n", seq_file.c_str());
            exit(-1);
        }

        init_fold(seq.c_str(), params);
        g_dangles = convert_to_dangle_mode(dangle_model);

        readThermodynamicParameters(param_dir.c_str());

        mpq_class energy = calculate(seq.length());

        write_header_subopt_file(output_file, seq, energy);

        ss_map_t subopt_data = subopt_traceback(seq.length(), delta, output_file, max_structure_count);

        free_fold(seq.length());
    };

    void generate_subopt(const std::string seq_file, const mpq_class delta, const ParameterVector params, const std::string param_dir, const int max_structure_count, const int dangle_model) {
        std::string output_file = subopt_filename(seq_file);
        generate_subopt(seq_file, output_file, delta, params, param_dir, max_structure_count, dangle_model);
    };
}
