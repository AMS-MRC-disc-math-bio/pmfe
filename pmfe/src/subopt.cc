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

namespace pmfe {
    using namespace std;

    static void write_header_subopt_file(std::string subopt_file, const std::string& seq, mpq_class energy)
    {
	ofstream outfile;
	outfile.open(subopt_file.c_str());
	char buff[4096];

	sprintf(buff,"count\t%s\t%f", seq.c_str(), energy.get_d());
	outfile << buff << std::endl;

	outfile.close();
    }

    static std::string subopt_filename(std::string seq_file, std::string output_dir) {
        std::string subopt_file = "";

        std::string outputPrefix = seq_file;

        int pos;
        // extract file name from the path
        if ((pos=outputPrefix.find_last_of('/')) > 0) {
            outputPrefix = outputPrefix.substr(pos+1);
        }

        // and if an extension exists, remove it ...
        if(outputPrefix.find(".") != std::string::npos)
            outputPrefix.erase(outputPrefix.rfind("."));

        // If output dir specified
        if (!output_dir.empty()) {
            subopt_file += output_dir;
            subopt_file += "/";
        }
        // ... and append the .ct
        subopt_file += outputPrefix;
        subopt_file += "_ss.txt";

        return subopt_file;
    };

    void generate_subopt(std::string seq_file, std::string output_dir, mpq_class delta, ParameterVector params, std::string param_dir, int max_structure_count, int dangle_model) {
        std::string seq = "";
        if (read_sequence_file(seq_file.c_str(), seq) == FAILURE) {
            printf("Failed to open sequence file: %s.\n\n", seq_file.c_str());
            exit(-1);
        }

        init_fold(seq.c_str(), params);
        if (dangle_model == 0 || dangle_model == 2)
            g_dangles = dangle_model;
        else
            exit(EXIT_FAILURE);

        readThermodynamicParameters(param_dir.c_str());

        std::string subopt_file = subopt_filename(seq_file, output_dir);

        mpq_class energy = calculate(seq.length()) ;
        write_header_subopt_file(subopt_file, seq, energy);

        ss_map_t subopt_data = subopt_traceback(seq.length(), delta, subopt_file, max_structure_count);

	free_fold(seq.length());
    };
}
