// Copyright (c) 2015 Andrew Gainer-Dewar

#include "catch.hpp"
#include <gmpxx.h>
#include <boost/filesystem.hpp>

#include "mfe.h"
#include "nndb_constants.h"
#include "nntm.h"
#include "pmfe_types.h"

namespace fs = boost::filesystem;

TEST_CASE("Combinatorial test sequence", "[mfe][combinatorial]") {
    // Load the sequence
    fs::path seqfile("test_data/combinatorial_seq.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 60);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants;
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE);
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-149, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.string() == "((((>..<((((....))))>..<((((....((((....))))....))))>..<))))");
    }

    SECTION("No multibranch penalties, positive d") {
        pmfe::ParameterVector params(0, 0, 0, 1);
        pmfe::Turner99 constants(params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE);
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-172, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.string() == "((((>..<((((....))))>..<((((....((((....))))....))))>..<))))");
    }

    SECTION("No multibranch penalties, negative d") {
        pmfe::ParameterVector params(0, 0, 0, -1);
        pmfe::Turner99 constants(params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE);
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-207, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.string() == "(..(....).)(....)(.(....)..)....(..(....).)(....)(.(....)..)");
    }

    SECTION("Artificial parameters") {
        pmfe::ParameterVector params(1, -1, 2, -2);
        pmfe::Turner99 constants(params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE);
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-414, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.string() == "(..(....).)(....)(.(....)..)....(..(....).)(....)(.(....)..)");
    }
}

TEST_CASE("Test tRNA sequence", "[mfe][test_tRNA]") {
    // Load the sequence
    fs::path seqfile("test_data/test_tRNA.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 76);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants;
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE);
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-229, 10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.string() == "(((((>....<((((((.((......)).))))))>.<((((..((((........))))..))))>)))))>...");
    }

    SECTION("No multibranch penalties, positive d") {
        pmfe::ParameterVector params(0, 0, 0, 1);
        pmfe::Turner99 constants(params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE);
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-139, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.string() == "((((((<(>((((........))))>((((((....).)))))><)><(((((.......)))))>))))))>...");
    }

    SECTION("No multibranch penalties, negative d") {
        pmfe::ParameterVector params(0, 0, 0, -1);
        pmfe::Turner99 constants(params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE);
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-432, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.string() == "(.(...)).((...)((...).)(...)(...)(..(.(.(...))))(...)(...)(....)(...))(...).");
    }

    SECTION("Artificial parameters") {
        pmfe::ParameterVector params(1, -1, 2, -2);
        pmfe::Turner99 constants(params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE);
        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);
        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-864, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.string() == "(.(...))..(...)((...).)(...)(...)(..(.(.(...))))(...)(...)(....)(...).(...).");
    }
}
