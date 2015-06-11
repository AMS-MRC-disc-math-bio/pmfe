/*
  GTfold: compute minimum free energy of RNA secondary structure
  Copyright (C) 2008  David A. Bader
  http://www.cc.gatech.edu/~bader

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <string.h>
#include <stdlib.h>
#include <utility>
#include <vector>
#include <iostream>
#include <stdexcept>

#include "nndb_constants.h"
#include "pmfe_types.h"
#include "thread_pool.h"
#include "rational.h"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/bind.hpp>

namespace pmfe
{

    namespace fs = boost::filesystem;

    std::vector< std::pair<RNA_base, RNA_base> > valid_pairs_in_order = {
        {BASE_A, BASE_U},
        {BASE_C, BASE_G},
        {BASE_G, BASE_C},
        {BASE_U, BASE_A},
        {BASE_G, BASE_U},
        {BASE_U, BASE_G}
    };

    std::vector< RNA_base > bases_in_order = {BASE_A, BASE_C, BASE_G, BASE_U};

    Turner99::Turner99(SimpleThreadPool& thread_pool, const ParameterVector& params, const fs::path& paramDir):
        NNDBConstants(params),
        thread_pool(thread_pool)
    {
        SimpleJobGroup job_group(thread_pool);

        job_group.post(boost::bind(&Turner99::initMiscValues, this, paramDir));

        job_group.post(boost::bind(&Turner99::initLoopValues, this, paramDir));

        job_group.post(boost::bind(&Turner99::initTstkhValues, this, paramDir));
        job_group.post(boost::bind(&Turner99::initTstkiValues, this, paramDir));

        job_group.post(boost::bind(&Turner99::initStackValues, this, paramDir));
        job_group.post(boost::bind(&Turner99::initDangleValues, this, paramDir));
        job_group.post(boost::bind(&Turner99::initTloopValues, this, paramDir));

        job_group.post(boost::bind(&Turner99::initIloop11Values, this, paramDir));
        job_group.post(boost::bind(&Turner99::initIloop21Values, this, paramDir));
        job_group.post(boost::bind(&Turner99::initIloop22Values, this, paramDir));

        job_group.wait_for_all_jobs();
    }

    void Turner99::initMiscValues(const fs::path& paramDir) {
        // Miscellaneous parameters
        fs::ifstream fileStream;

        fs::path fileName = paramDir / "miscloop.DAT";
        fileStream.open(fileName);

        if (!fileStream.good()) {
            std::stringstream errorMessage;
            errorMessage << "Parameter file " << fileName << "is invalid.";
            throw std::runtime_error(errorMessage.str());
        }

        std::string currentWord;

        // Extrapolation for large loops based on polymer theory
        // internal, bulge or hairpin loops > 30: dS(T) = dS(30) + prelog * ln(n/30)
        fileStream >> currentWord;
        prelog = params.dummy_scaling * get_rational_from_word(currentWord);

        // Maximum correction for asymmetric internal loops per Ninio
        fileStream >> currentWord;
        maxpen = params.dummy_scaling * get_rational_from_word(currentWord);

        // Ninio's f(m) array
        for (int count = 1; count <= 4; ++count) {
            fileStream >> currentWord;
            poppen[count] = params.dummy_scaling * get_rational_from_word(currentWord);
        }

        // Multibranch loop
        // Skip the lines in the file and set these values manually
        fileStream >> currentWord;
        multConst[0] = params.multiloop_penalty;

        fileStream >> currentWord;
        multConst[1] = params.unpaired_penalty;

        fileStream >> currentWord;
        multConst[2] = params.branch_penalty;

        // efn values, unused in our model
        fileStream >> currentWord;
        fileStream >> currentWord;
        fileStream >> currentWord;

        // Non-GC terminal pair penalty (sometimes called "AU penalty")
        fileStream >> currentWord;
        auend = params.dummy_scaling * get_rational_from_word(currentWord);

        // Bonus for GGG hairpin
        fileStream >> currentWord;
        gubonus = params.dummy_scaling * get_rational_from_word(currentWord);

        // C hairpin slope
        fileStream >> currentWord;
        cslope = params.dummy_scaling * (get_rational_from_word(currentWord) + Rational(1,100));

        // C hairpin intercept
        fileStream >> currentWord;
        cint = params.dummy_scaling * get_rational_from_word(currentWord);

        // C hairpin of 3
        fileStream >> currentWord;
        c3 = params.dummy_scaling * (get_rational_from_word(currentWord) + Rational(1, 100));

        // init value, unused in our model
        fileStream >> currentWord;

        // GAIL rule (Grossly Asymmetric Interior Loop)
        fileStream >> currentWord;
        gail = (currentWord != "0");
    }

    void Turner99::initLoopValues(const fs::path& paramDir) {
        // Destabilizing energies of various types of loops
        fs::ifstream fileStream;

        fs::path fileName = paramDir / "loop.DAT";
        fileStream.open(fileName);

        if (!fileStream.good()) {
            std::stringstream errorMessage;
            errorMessage << "Parameter file " << fileName << "is invalid.";
            throw std::runtime_error(errorMessage.str());
        }

        for (int loopsize = 1; loopsize <= 30; ++loopsize) {
            // Each size of loop gets its own line in the file
            std::string line;
            while (line.empty()) {
                getline(fileStream, line);
            }

            std::stringstream linestream;
            linestream << line;
            std::string currentWord;

            // Column 1: the index, which we already know
            linestream >> currentWord;

            // Column 2: internal loop
            linestream >> currentWord;
            if (currentWord == "inf") {
                inter[loopsize] = Rational::infinity();
            } else {
                inter[loopsize] = params.dummy_scaling * get_rational_from_word(currentWord);
            }

            // Column 3: bulge loop
            linestream >> currentWord;
            if (currentWord == "inf") {
                bulge[loopsize] = Rational::infinity();
            } else {
                bulge[loopsize] = params.dummy_scaling * get_rational_from_word(currentWord);
            }

            // Column 4: hairpin loop
            linestream >> currentWord;
            if (currentWord == "inf") {
                hairpin[loopsize] = Rational::infinity();
            } else {
                hairpin[loopsize] = params.dummy_scaling * get_rational_from_word(currentWord);
            }
        }
    }

    void Turner99::initTstkhValues(const fs::path& paramDir) {
        /*
          tstkh: Terminal mismatches in hairpin loops
          tstkh[A][B][C][D] represents the energy associated with this configuration:
          5' ---> 3'
              AC
              BD
          3' <--- 5'
        */

        fs::ifstream fileStream;

        fs::path fileName = paramDir / "tstackh.DAT";
        fileStream.open(fileName);

        if (!fileStream.good()) {
            std::stringstream errorMessage;
            errorMessage << "Parameter file " << fileName << "is invalid.";
            throw std::runtime_error(errorMessage.str());
        }

        for (const auto& A: bases_in_order) {
            for (const auto& C: bases_in_order) {
                // Get the line of parameters for these bases A and C
                std::string line;
                while (line.empty()) {
                    getline(fileStream, line);
                }

                std::stringstream linestream;
                linestream << line;
                std::string currentWord;

                for (const auto& B: bases_in_order) {
                    for (const auto& D: bases_in_order) {
                        linestream >> currentWord;
                        if (currentWord == "inf") {
                            tstkh[A][B][C][D] = Rational::infinity();
                        } else {
                            tstkh[A][B][C][D] = params.dummy_scaling * get_rational_from_word(currentWord);
                        }
                    }
                }
            }
        }
    }

    void Turner99::initTstkiValues(const fs::path& paramDir) {
        /*
          tstki: Terminal mismatches in internal loops
          tstki[A][B][C][D] represents the energy associated with this configuration:
          5' ---> 3'
              AC
              BD
          3' <--- 5'
        */

        fs::ifstream fileStream;

        fs::path fileName = paramDir / "tstacki.DAT";
        fileStream.open(fileName);

        if (!fileStream.good()) {
            std::stringstream errorMessage;
            errorMessage << "Parameter file " << fileName << "is invalid.";
            throw std::runtime_error(errorMessage.str());
        }

        for (const auto& A: bases_in_order) {
            for (const auto& C: bases_in_order) {
                // Get the line of parameters for these bases A and C
                std::string line;
                while (line.empty()) {
                    getline(fileStream, line);
                }

                std::stringstream linestream;
                linestream << line;
                std::string currentWord;

                for (const auto& B: bases_in_order) {
                    for (const auto& D: bases_in_order) {
                        linestream >> currentWord;
                        if (currentWord == "inf") {
                            tstki[A][B][C][D] = Rational::infinity();
                        } else {
                            tstki[A][B][C][D] = params.dummy_scaling * get_rational_from_word(currentWord);
                        }
                    }
                }
            }
        }
    }

    void Turner99::initStackValues(const fs::path& paramDir) {
        /*
          stack: Energies associated with stacks
          stack[A][B][C][D] represents the energy associated with this configuration:
          5' ---> 3'
              AC
              BD
          3' <--- 5'
        */

        fs::ifstream fileStream;

        fs::path fileName = paramDir / "stack.DAT";
        fileStream.open(fileName);

        if (!fileStream.good()) {
            std::stringstream errorMessage;
            errorMessage << "Parameter file " << fileName << "is invalid.";
            throw std::runtime_error(errorMessage.str());
        }

        for (const auto& A: bases_in_order) {
            for (const auto& C: bases_in_order) {
                // Get the line of parameters for these bases A and C
                std::string line;
                while (line.empty()) {
                    getline(fileStream, line);
                }

                std::stringstream linestream;
                linestream << line;
                std::string currentWord;

                for (const auto& B: bases_in_order) {
                    for (const auto& D: bases_in_order) {
                        linestream >> currentWord;
                        if (currentWord == "inf") {
                            stack[A][B][C][D] = Rational::infinity();
                        } else {
                            stack[A][B][C][D] = params.dummy_scaling * get_rational_from_word(currentWord);
                        }
                    }
                }
            }
        }
    }

    void Turner99::initDangleValues(const fs::path& paramDir) {
        /*
          dangle: Energies associated with dangling bases
          dangle[A][B][C][0] represents the energy associated with this configuration, where
          C is a "3' dangling end":
          5' ---> 3'
              AC
              B
          3' <--- 5'

          dangle[A][B][C][1] represents the energy associated with this configuration, where
          C is a "5' dangling end":
          5' ---> 3'
              A
              BC
          3' <--- 5'
        */

        fs::ifstream fileStream;

        fs::path fileName = paramDir / "dangle.DAT";
        fileStream.open(fileName);

        if (!fileStream.good()) {
            std::stringstream errorMessage;
            errorMessage << "Parameter file " << fileName << "is invalid.";
            throw std::runtime_error(errorMessage.str());
        }

        for (int m = 0; m <= 1; ++m) {
            for (const auto& A: bases_in_order) {
                // Get the line of parameters for this base A and mode m
                std::string line;
                while (line.empty()) {
                    getline(fileStream, line);
                }

                std::stringstream linestream;
                linestream << line;
                std::string currentWord;

                for (const auto& B: bases_in_order) {
                    for (const auto& C: bases_in_order) {
                        linestream >> currentWord;
                        if (currentWord == "inf") {
                            dangle[A][B][C][m] = Rational::infinity();
                        } else {
                            dangle[A][B][C][m] = params.dummy_scaling * get_rational_from_word(currentWord);
                        }
                    }
                }
            }
        }
    }

    void Turner99::initTloopValues(const fs::path& paramDir) {
        // Custom energies for certain short hairpins
        fs::ifstream fileStream;

        fs::path fileName = paramDir / "tloop.DAT";
        fileStream.open(fileName);

        if (!fileStream.good()) {
            std::stringstream errorMessage;
            errorMessage << "Parameter file " << fileName << "is invalid.";
            throw std::runtime_error(errorMessage.str());
        }

        while (fileStream.good() and not fileStream.eof()) {
            // Each custom loop has its own line in the file
            std::string line;
            getline(fileStream, line);
            if (line.empty()) {
                continue;
            }

            std::stringstream linestream;
            linestream << line;
            std::string currentWord;

            // Loop as a subsequence
            linestream >> currentWord;
            std::string loopkey = currentWord;

            // Associated energy
            linestream >> currentWord;
            Rational loopval = params.dummy_scaling * get_rational_from_word(currentWord);

            tloop[loopkey] = loopval;
        }
    }

    void Turner99::initIloop11Values(const fs::path& paramDir) {
        /*
          iloop11: Energies associated with internal loops with one base on each side

          iloop11[A][B][C][D][E][F] represents the energy associated with this
          configuration, where AD and FILESTREAM are pairs and B and E form the loop:

          5' ----> 3'
               B
              A C
              D F
               E
          3' <---- 5'
        */

        fs::ifstream fileStream;

        fs::path fileName = paramDir / "int11.DAT";
        fileStream.open(fileName);

        if (!fileStream.good()) {
            std::stringstream errorMessage;
            errorMessage << "Parameter file " << fileName << "is invalid.";
            throw std::runtime_error(errorMessage.str());
        }

        // Pre-populate the array with Rational::infinity()
        for (const auto& A: bases_in_order)
            for (const auto& B: bases_in_order)
                for (const auto& C: bases_in_order)
                    for (const auto& D: bases_in_order)
                        for (const auto& E: bases_in_order)
                            for (const auto& F: bases_in_order)
                                iloop11[A][B][C][D][E][F] = Rational::infinity();

        // Then read the parameters from the file
        for (const auto& left_pair: valid_pairs_in_order) {
            RNA_base A = left_pair.first;
            RNA_base D = left_pair.second;
            for (const RNA_base& B: bases_in_order) {
                // Read the line of parameters for these bases A, B, and D
                std::string line;
                while (line.empty()) {
                    getline(fileStream, line);
                }

                std::stringstream linestream;
                linestream << line;
                std::string currentWord;

                for (const auto& right_pair: valid_pairs_in_order) {
                    RNA_base C = right_pair.first;
                    RNA_base F = right_pair.second;
                    for (const RNA_base& E: bases_in_order) {
                        linestream >> currentWord;
                        iloop11[A][B][C][D][E][F] = params.dummy_scaling * get_rational_from_word(currentWord);
                    }
                }
            }
        }
    }

    void Turner99::initIloop21Values(const fs::path& paramDir) {
        /*
          iloop21: Energies associated with internal loops with one base
          on one side and two on the other

          iloop21[A][B][C][D][E][F][G] represents the energy associated with this
          configuration, where AB and FG are pairs, C is the 1-base side of the loop,
          and DE is the 2-base side:

          5' ----> 3'
              C
             A  F
             B  G
              DE
          3' <---- 5'
        */

        fs::ifstream fileStream;

        fs::path fileName = paramDir / "int21.DAT";
        fileStream.open(fileName);

        if (!fileStream.good()) {
            std::stringstream errorMessage;
            errorMessage << "Parameter file " << fileName << "is invalid.";
            throw std::runtime_error(errorMessage.str());
        }

        // Pre-populate the array with Rational::infinity()
        for (const auto& A: bases_in_order)
            for (const auto& B: bases_in_order)
                for (const auto& C: bases_in_order)
                    for (const auto& D: bases_in_order)
                        for (const auto& E: bases_in_order)
                            for (const auto& F: bases_in_order)
                                for (const auto& G: bases_in_order)
                                    iloop21[A][B][C][D][E][F][G] = Rational::infinity();

        // Then read the parameters from the file
        for (const auto& left_pair: valid_pairs_in_order) {
            RNA_base A = left_pair.first;
            RNA_base B = left_pair.second;
            for (const RNA_base& E: bases_in_order) {
                for (const RNA_base& C: bases_in_order) {
                    // Read the line of parameters for these bases A, B, C, and E
                    std::string line;
                    while (line.empty()) {
                        getline(fileStream, line);
                    }

                    std::stringstream linestream;
                    linestream << line;
                    std::string currentWord;

                    for (const auto& right_pair: valid_pairs_in_order) {
                        RNA_base F = right_pair.first;
                        RNA_base G = right_pair.second;
                        for (const RNA_base& D: bases_in_order) {
                            linestream >> currentWord;
                            iloop21[A][B][C][D][E][F][G] = params.dummy_scaling * get_rational_from_word(currentWord);
                        }
                    }
                }
            }
        }
    }


    void Turner99::initIloop22Values(const fs::path& paramDir) {
        /*
          iloop22: Energies associated with internal loops with two bases on each side

          iloop22[A][B][C][D][E][F][G][H] represents the energy associated with this
          configuration, where AC and BD are pairs and E, F, G, and H form the loop:
          5' ----> 3'
              EG
             A  B
             C  D
              FH
          3' <---- 5'
        */

        fs::ifstream fileStream;

        fs::path fileName = paramDir / "int22.DAT";
        fileStream.open(fileName);

        if (!fileStream.good()) {
            std::stringstream errorMessage;
            errorMessage << "Parameter file " << fileName << "is invalid.";
            throw std::runtime_error(errorMessage.str());
        }

        // Pre-populate the array with Rational::infinity()
        for (const auto& A: bases_in_order)
            for (const auto& B: bases_in_order)
                for (const auto& C: bases_in_order)
                    for (const auto& D: bases_in_order)
                        for (const auto& E: bases_in_order)
                            for (const auto& F: bases_in_order)
                                for (const auto& G: bases_in_order)
                                    for (const auto& H: bases_in_order)
                                        iloop22[A][B][C][D][E][F][G][H] = Rational::infinity();

        // Then read the parameters from the file
        for (const auto& left_pair: valid_pairs_in_order) {
            RNA_base A = left_pair.first;
            RNA_base C = left_pair.second;
            for (const auto& right_pair: valid_pairs_in_order) {
                RNA_base B = right_pair.first;
                RNA_base D = right_pair.second;
                for (const RNA_base& E: bases_in_order) {
                    for (const RNA_base& F: bases_in_order) {
                        // Read the line of parameters
                        std::string line;
                        while (line.empty()) {
                            getline(fileStream, line);
                        }

                        std::stringstream linestream;
                        linestream << line;
                        std::string currentWord;

                        for (const RNA_base& G: bases_in_order) {
                            for (const RNA_base& H: bases_in_order) {
                                linestream >> currentWord;
                                iloop22[A][B][C][D][E][F][G][H] = params.dummy_scaling * get_rational_from_word(currentWord);
                            }
                        }
                    }
                }
            }
        }
    }
}
