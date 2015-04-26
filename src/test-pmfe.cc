// Copyright (c) 2015 Andrew Gainer-Dewar

#include "catch.hpp"
#include <gmpxx.h>
#include <boost/filesystem.hpp>

#include "mfe.h"
#include "nndb_constants.h"
#include "nntm.h"
#include "pmfe_types.h"

namespace fs = boost::filesystem;

TEST_CASE("Combinatorial test sequence MFE", "[mfe][combinatorial]") {
    // Load the sequence
    fs::path seqfile("test_data/combinatorial_seq.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 60);

    SECTION("Turner99 published parameters, d1") {
        pmfe::Turner99 constants;
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE);
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-149, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.string() == "((((>..<((((....))))>..<((((....((((....))))....))))>..<))))");

        mpq_class delta = 1;

        std::vector<pmfe::RNAStructureWithScore> suboptimal_structures = energy_model.suboptimal_structures(seq_annotated, delta, true);

        REQUIRE(suboptimal_structures.size() == 34);
        REQUIRE(suboptimal_structures.front().score.energy == energy);
        REQUIRE(suboptimal_structures.back().score.energy <= energy + delta);
    }

    SECTION("Turner99 published parameters, d2") {
        pmfe::Turner99 constants;
        pmfe::NNTM energy_model(constants, pmfe::BOTH_DANGLE);
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-149, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.string() == "((((....((((....))))....((((....((((....))))....))))....))))");

        mpq_class delta = 1;

        std::vector<pmfe::RNAStructureWithScore> suboptimal_structures = energy_model.suboptimal_structures(seq_annotated, delta, true);

        REQUIRE(suboptimal_structures.size() == 5);
        REQUIRE(suboptimal_structures.front().score.energy == energy);
        REQUIRE(suboptimal_structures.back().score.energy <= energy + delta);
    }

    SECTION("No multibranch penalties, positive d, d1") {
        pmfe::ParameterVector params(0, 0, 0, 1);
        pmfe::Turner99 constants(params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE);
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-172, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.string() == "((((>..<((((....))))>..<((((....((((....))))....))))>..<))))");

        mpq_class delta = 1;

        std::vector<pmfe::RNAStructureWithScore> suboptimal_structures = energy_model.suboptimal_structures(seq_annotated, delta, true);

        REQUIRE(suboptimal_structures.size() == 39);
        REQUIRE(suboptimal_structures.front().score.energy == energy);
        REQUIRE(suboptimal_structures.back().score.energy <= energy + delta);
    }

    SECTION("No multibranch penalties, positive d, d2") {
        pmfe::ParameterVector params(0, 0, 0, 1);
        pmfe::Turner99 constants(params);
        pmfe::NNTM energy_model(constants, pmfe::BOTH_DANGLE);
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-172, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.string() == "((((....((((....))))....((((....((((....))))....))))....))))");

        mpq_class delta = 1;

        std::vector<pmfe::RNAStructureWithScore> suboptimal_structures = energy_model.suboptimal_structures(seq_annotated, delta, true);

        REQUIRE(suboptimal_structures.size() == 6);
        REQUIRE(suboptimal_structures.front().score.energy == energy);
        REQUIRE(suboptimal_structures.back().score.energy <= energy + delta);
    }

    SECTION("No multibranch penalties, negative d, d1") {
        pmfe::ParameterVector params(0, 0, 0, -1);
        pmfe::Turner99 constants(params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE);
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-207, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.string() == "(..(....).)(....)(.(....)..)....(..(....).)(....)(.(....)..)");

        mpq_class delta = 2;

        std::vector<pmfe::RNAStructureWithScore> suboptimal_structures = energy_model.suboptimal_structures(seq_annotated, delta, true);

        REQUIRE(suboptimal_structures.size() == 82);
        REQUIRE(suboptimal_structures.front().score.energy == energy);
        REQUIRE(suboptimal_structures.back().score.energy <= energy + delta);
    }

    SECTION("No multibranch penalties, negative d, d2") {
        pmfe::ParameterVector params(0, 0, 0, -1);
        pmfe::Turner99 constants(params);
        pmfe::NNTM energy_model(constants, pmfe::BOTH_DANGLE);
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-35));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.string() == "(..(....).)(....)(.(....)..)....(..(....).)(....)(.(....)..)");

        mpq_class delta = 2;

        std::vector<pmfe::RNAStructureWithScore> suboptimal_structures = energy_model.suboptimal_structures(seq_annotated, delta, true);

        REQUIRE(suboptimal_structures.size() == 24);
        REQUIRE(suboptimal_structures.front().score.energy == energy);
        REQUIRE(suboptimal_structures.back().score.energy <= energy + delta);
   }

    SECTION("Artificial parameters, d1") {
        pmfe::ParameterVector params(1, -1, 2, -2);
        pmfe::Turner99 constants(params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE);
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-414, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.string() == "(..(....).)(....)(.(....)..)....(..(....).)(....)(.(....)..)");

        mpq_class delta = 2;

        std::vector<pmfe::RNAStructureWithScore> suboptimal_structures = energy_model.suboptimal_structures(seq_annotated, delta, true);

        REQUIRE(suboptimal_structures.size() == 4);
        REQUIRE(suboptimal_structures.front().score.energy == energy);
        REQUIRE(suboptimal_structures.back().score.energy <= energy + delta);
    }

    SECTION("Artificial parameters, d2") {
        pmfe::ParameterVector params(1, -1, 2, -2);
        pmfe::Turner99 constants(params);
        pmfe::NNTM energy_model(constants, pmfe::BOTH_DANGLE);
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-70));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.string() == "(..(....).)(....)(.(....)..)....(..(....).)(....)(.(....)..)");

        mpq_class delta = 2;

        std::vector<pmfe::RNAStructureWithScore> suboptimal_structures = energy_model.suboptimal_structures(seq_annotated, delta, true);

        REQUIRE(suboptimal_structures.size() == 3);
        REQUIRE(suboptimal_structures.front().score.energy == energy);
        REQUIRE(suboptimal_structures.back().score.energy <= energy + delta);
    }
}

TEST_CASE("Test tRNA sequence MFE", "[mfe][test_tRNA]") {
    // Load the sequence
    fs::path seqfile("test_data/test_tRNA.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 76);

    SECTION("Turner99 published parameters, d1") {
        pmfe::Turner99 constants;
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE);
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-229, 10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.string() == "(((((>....<((((((.((......)).))))))>.<((((..((((........))))..))))>)))))>...");

        mpq_class delta = 1;

        std::vector<pmfe::RNAStructureWithScore> suboptimal_structures = energy_model.suboptimal_structures(seq_annotated, delta, true);

        REQUIRE(suboptimal_structures.size() == 55);
        REQUIRE(suboptimal_structures.front().score.energy == energy);
        REQUIRE(suboptimal_structures.back().score.energy <= energy + delta);
    }

    SECTION("Turner99 published parameters, d2") {
        pmfe::Turner99 constants;
        pmfe::NNTM energy_model(constants, pmfe::BOTH_DANGLE);
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-119, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.string() == "((((((.....((((((.((......)).))))))...((((..((((........))))..))))))))))....");

        mpq_class delta = 1;

        std::vector<pmfe::RNAStructureWithScore> suboptimal_structures = energy_model.suboptimal_structures(seq_annotated, delta, true);

        REQUIRE(suboptimal_structures.size() == 14);
        REQUIRE(suboptimal_structures.front().score.energy == energy);
        REQUIRE(suboptimal_structures.back().score.energy <= energy + delta);
    }

    SECTION("No multibranch penalties, positive d, d1") {
        pmfe::ParameterVector params(0, 0, 0, 1);
        pmfe::Turner99 constants(params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE);
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-139, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.string() == "((((((<(>((((........))))>((((((....).)))))><)><(((((.......)))))>))))))>...");

        mpq_class delta = 1;

        std::vector<pmfe::RNAStructureWithScore> suboptimal_structures = energy_model.suboptimal_structures(seq_annotated, delta, true);

        REQUIRE(suboptimal_structures.size() == 257);
        REQUIRE(suboptimal_structures.front().score.energy == energy);
        REQUIRE(suboptimal_structures.back().score.energy <= energy + delta);
    }

    SECTION("No multibranch penalties, positive d, d2") {
        pmfe::ParameterVector params(0, 0, 0, 1);
        pmfe::Turner99 constants(params);
        pmfe::NNTM energy_model(constants, pmfe::BOTH_DANGLE);
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-158, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.string() == "((((((((.((((........)))).((((((....).))))))))))(((((.......)))))....)))....");

        mpq_class delta = 2;

        std::vector<pmfe::RNAStructureWithScore> suboptimal_structures = energy_model.suboptimal_structures(seq_annotated, delta, true);

        REQUIRE(suboptimal_structures.size() == 26);
        REQUIRE(suboptimal_structures.front().score.energy == energy);
        REQUIRE(suboptimal_structures.back().score.energy <= energy + delta);

    }

    SECTION("No multibranch penalties, negative d, d1") {
        pmfe::ParameterVector params(0, 0, 0, -1);
        pmfe::Turner99 constants(params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE);
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-432, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.string() == "(.(...)).((...)((...).)(...)(...)(..(.(.(...))))(...)(...)(....)(...))(...).");

        mpq_class delta = 1;

        std::vector<pmfe::RNAStructureWithScore> suboptimal_structures = energy_model.suboptimal_structures(seq_annotated, delta, true);

        REQUIRE(suboptimal_structures.size() == 124);
        REQUIRE(suboptimal_structures.front().score.energy == energy);
        REQUIRE(suboptimal_structures.back().score.energy <= energy + delta);
    }

    SECTION("No multibranch penalties, negative d, d2") {
        pmfe::ParameterVector params(0, 0, 0, -1);
        pmfe::Turner99 constants(params);
        pmfe::NNTM energy_model(constants, pmfe::BOTH_DANGLE);
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-787, 10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.string() == "(.(...))(.(...))(.((..(..(..(.(.(((...).)..)..)..).).).)..)..).(...)(.(...))");

        mpq_class delta = 1;

        std::vector<pmfe::RNAStructureWithScore> suboptimal_structures = energy_model.suboptimal_structures(seq_annotated, delta, true);

        REQUIRE(suboptimal_structures.size() == 92);
        REQUIRE(suboptimal_structures.front().score.energy == energy);
        REQUIRE(suboptimal_structures.back().score.energy <= energy + delta);
    }

    SECTION("Artificial parameters, d1") {
        pmfe::ParameterVector params(1, -1, 2, -2);
        pmfe::Turner99 constants(params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE);
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-864, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.string() == "(.(...))..(...)((...).)(...)(...)(..(.(.(...))))(...)(...)(....)(...).(...).");

        mpq_class delta = 2;

        std::vector<pmfe::RNAStructureWithScore> suboptimal_structures = energy_model.suboptimal_structures(seq_annotated, delta, true);

        REQUIRE(suboptimal_structures.size() == 53);
        REQUIRE(suboptimal_structures.front().score.energy == energy);
        REQUIRE(suboptimal_structures.back().score.energy <= energy + delta);
    }

    SECTION("Artificial parameters, d2") {
        pmfe::ParameterVector params(1, -1, 2, -2);
        pmfe::Turner99 constants(params);
        pmfe::NNTM energy_model(constants, pmfe::BOTH_DANGLE);
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-787, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.string() == "(.(...))(.(...))(.((..(..(..(.(.(((...).)..)..)..).).).)..)..).(...)(.(...))");

        mpq_class delta = 2;

        std::vector<pmfe::RNAStructureWithScore> suboptimal_structures = energy_model.suboptimal_structures(seq_annotated, delta, true);

        REQUIRE(suboptimal_structures.size() == 82);
        REQUIRE(suboptimal_structures.front().score.energy == energy);
        REQUIRE(suboptimal_structures.back().score.energy <= energy + delta);
    }
}
