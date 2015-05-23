// Copyright (c) 2015 Andrew Gainer-Dewar

#include "catch.hpp"
#include <gmpxx.h>
#include <boost/filesystem.hpp>

#include "mfe.h"
#include "nndb_constants.h"
#include "nntm.h"
#include "pmfe_types.h"
#include "thread_pool.h"

namespace fs = boost::filesystem;

TEST_CASE("A. tabira 5S MFE", "[mfe][biological][atabira][5S]") {
    // Build the thread pool
    pmfe::SimpleThreadPool thread_pool;

    // Load the sequence
    fs::path seqfile("test_seq/5S/a.tabira_5S.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 120);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-533, 10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((((.....((((((((((...((.((((......))))..))..)))...)))))))..(((((((...(.((..(..((....))..)..)))..))))))))))))))))..");
    }
}

TEST_CASE("C. diphtheriae tRNA MFE", "[mfe][biological][cdiphtheriae][tRNA]") {
    // Build the thread pool
    pmfe::SimpleThreadPool thread_pool;

    // Load the sequence
    fs::path seqfile("test_seq/tRNA/c.diphtheriae_tRNA.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 74);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-35));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "((((((((((.(.(.((((...((.(((.(((((((....)))).))).))).)))))).))))).))))))).");
    }
}

TEST_CASE("D. mobilis 5S MFE", "[mfe][biological][dmobilis][5S]") {
    // Build the thread pool
    pmfe::SimpleThreadPool thread_pool;

    // Load the sequence
    fs::path seqfile("test_seq/5S/d.mobilis_5S.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 133);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-463, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "..(((((((((((((((....((((((((..(((..((((..(((...)))..))))..)))...)))))))).((((((((((.((((((.((....)))))))).)))))))))))))))).)))))))))");
    }
}

TEST_CASE("E. coli 5S MFE", "[mfe][biological][ecoli][5S]") {
    // Build the thread pool
    pmfe::SimpleThreadPool thread_pool;

    // Load the sequence
    fs::path seqfile("test_seq/5S/e.coli_5S.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 120);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-549,10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "((((((((((.........((((....)))).((((((((((.....(((....)))...(((((.......))))).))))))))))..(((.(((....))))))..)))))))))).");
    }
}

TEST_CASE("G. arboreum 5S MFE", "[mfe][biological][garboreum][5S]") {
    // Build the thread pool
    pmfe::SimpleThreadPool thread_pool;

    // Load the sequence
    fs::path seqfile("test_seq/5S/g.arboreum_5S.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 120);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-407, 10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((((...((.((.(((....................(((((..(.((((((....)))))).).)))))...((((((.((....))))))))..))).)).)))))))))))..");
    }
}

TEST_CASE("H. sapiens tRNA MFE", "[mfe][biological][hsapiens][tRNA]") {
    // Build the thread pool
    pmfe::SimpleThreadPool thread_pool;

    // Load the sequence
    fs::path seqfile("test_seq/tRNA/h.sapiens_tRNA.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 72);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-263, 10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "((((.....)))).((((((.((.(((((((((..((((....))))..))).))))))))...))))))..");
    }
}

TEST_CASE("L. delbrueckii tRNA MFE", "[mfe][biological][ldelbrueckii][tRNA]") {
    // Build the thread pool
    pmfe::SimpleThreadPool thread_pool;

    // Load the sequence
    fs::path seqfile("test_seq/tRNA/l.delbrueckii_tRNA.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 72);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-239, 10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((...............(((.(((((.......))))).)))..(((.......)))..))))))).");
    }
}

TEST_CASE("O. nivara tRNA MFE", "[mfe][biological][onivara][tRNA]") {
    // Build the thread pool
    pmfe::SimpleThreadPool thread_pool;

    // Load the sequence
    fs::path seqfile("test_seq/tRNA/o.nivara_tRNA.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 73);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-271, 10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((..((((........))))((((((.......))))))....(((((.......)))))))))))).");
    }
}

TEST_CASE("R. norvegicus 5S MFE", "[mfe][biological][rnorvegicus][5S]") {
    // Build the thread pool
    pmfe::SimpleThreadPool thread_pool;

    // Load the sequence
    fs::path seqfile("test_seq/5S/r.norvegicus_5S.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 123);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-539, 10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((((....((((((((...((..((((..((....))..))))..))....)))))).))(((((((...(.((..(.(((....))).)..)))..)))))))))))))))).....");
    }
}

TEST_CASE("S. tokodaii tRNA MFE", "[mfe][biological][stokodaii][tRNA]") {
    // Build the thread pool
    pmfe::SimpleThreadPool thread_pool;

    // Load the sequence
    fs::path seqfile("test_seq/tRNA/s.tokodaii_tRNA.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 74);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-79, 2));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "((((((((.((((.((...)).))))..((((.(((.((((((.......)))))).))).)))))))))))).");
    }
}

TEST_CASE("Combinatorial sequence MFE", "[mfe][synthetic][combinatorial]") {
    // Build the thread pool
    pmfe::SimpleThreadPool thread_pool;

    // Load the sequence
    fs::path seqfile("test_seq/synthetic/test_combinatorial.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 60);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-149, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "((((....((((....))))....((((....((((....))))....))))....))))");
    }
}

TEST_CASE("Randomly generated sequence MFE", "[mfe][synthetic][random]") {
    // Build the thread pool
    pmfe::SimpleThreadPool thread_pool;

    // Load the sequence
    fs::path seqfile("test_seq/synthetic/test_random.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 75);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-81, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                ".(((.(((....(((....((.....))..))).....))))))................((......)).....");
    }
}

TEST_CASE("O. nivara tRNA (old) MFE", "[mfe][biological][onivara-old][trna]") {
    // Build the thread pool
    pmfe::SimpleThreadPool thread_pool;

    // Load the sequence
    std::string theseq = "AUCAGAGUGGCGCAGCGGAAGCGUGGUGGGCCCAUAACCCACAGGUCCCAGGAUCGAAACCUGGCUCUGAUA";
    pmfe::RNASequence seq(theseq);

    // Some basic sanity checks
    REQUIRE(seq.len() == 72);

    SECTION("onivara_old_classical, d1") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-283, 10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "((((((((..(((.((....)))))(((((.......))))).((((....))))........)))))))).");
    }

    SECTION("onivara_old_1, d1") {
        pmfe::ParameterVector params(1, 1, 1, 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-239, 10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((.(((.((((....)).))(((((.......)))))..)))(((((.......)))))))))))).");
    }

    SECTION("onivara_old_2, d1") {
        pmfe::ParameterVector params(-1, 1, 1, 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-697, 25));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((.(((.((((....)).))(((((.......)))))..)))(((((.......)))))))))))).");
    }

    SECTION("onivara_old_3, d1") {
        pmfe::ParameterVector params(-1, 0, 1, 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-1669, 50));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((.(((...((....))...(((((.......)))))..)))(((((.......)))))))))))).");
    }

    SECTION("onivara_old_4, d1") {
        pmfe::ParameterVector params(1, 0, -1, 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-2067, 50));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((.(((...((....))...(((((.......)))))..)))(((((.......)))))))))))).");
    }

    SECTION("onivara_old_5, d1") {
        pmfe::ParameterVector params(mpq_class(1, 10), 1, mpq_class(-1, 10), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-806, 25));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((.(((.((((....)).))(((((.......)))))..)))(((((.......)))))))))))).");
    }

    SECTION("onivara_old_6, d1") {
        pmfe::ParameterVector params(3, mpq_class(-1, 10), 1, 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-2729, 100));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "((((((((..(((.((....)))))(((((.......))))).((((....))))........)))))))).");
    }

    SECTION("onivara_old_7, d1") {
        pmfe::ParameterVector params(2, mpq_class(1, 10), 5, 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-108, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "..........(((.((....)))))(((((.......))))).((..(((((.......)))))..))....");
    }

    SECTION("onivara_old_8, d1") {
        pmfe::ParameterVector params(2, mpq_class(1, 10), 0, 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-65, 2));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((.(((...((....))...(((((.......)))))..)))(((((.......)))))))))))).");
    }

    SECTION("onivara_old_9, d1") {
        pmfe::ParameterVector params(0, mpq_class(1, 10), 0, 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-73, 2));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((.(((...((....))...(((((.......)))))..)))(((((.......)))))))))))).");
    }

    SECTION("onivara_old_10, d1") {
        pmfe::ParameterVector params(mpq_class(-1, 10), mpq_class(1, 10), 2, 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-617, 25));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((.(((...((....))...(((((.......)))))..)))(((((.......)))))))))))).");
    }

    SECTION("onivara_old_11, d1") {
        pmfe::ParameterVector params(mpq_class(17, 5), mpq_class(1, 10), mpq_class(2, 5), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-55, 2));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((((((((.((....)))))(((((.......))))).((((....))))....))..)))))))).");
    }


    SECTION("onivara_old_12, d1") {
        pmfe::ParameterVector params(mpq_class(17, 5), mpq_class(1, 5), mpq_class(2, 5), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-27, 1));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((((((((.((....)))))(((((.......))))).((((....))))....))..)))))))).");
    }

    SECTION("onivara_old_13, d1") {
        pmfe::ParameterVector params(mpq_class(16, 5), mpq_class(1, 5), mpq_class(2, 5), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-136, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((((((((.((....)))))(((((.......))))).((((....))))....))..)))))))).");
    }

    SECTION("onivara_old_14, d1") {
        pmfe::ParameterVector params(mpq_class(16, 5), 0, mpq_class(2, 5), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-5034, 25));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "((....)((((...)(....))((...)((.......)((...)))))(((...)(....))(....)))).");
    }

    SECTION("onivara_old_15, d1") {
        pmfe::ParameterVector params(mpq_class(43, 5), mpq_class(31, 10), mpq_class(-131, 10), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-143, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((.(((...((....))...(((((.......)))))..)))(((((.......)))))))))))).");
    }

    SECTION("onivara_old_16, d1") {
        pmfe::ParameterVector params(mpq_class(18, 5), mpq_class(-63, 10), 2, 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-35461, 100));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(..........(...)..................................(...)...............).");
    }

    SECTION("onivara_17, d1") {
        pmfe::ParameterVector params(mpq_class(-7, 5), mpq_class(39, 10), -10, 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-5172, 25));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "((....)((((...)(....))((...)((.......)((...)))))(((...)(....))(....)))).");
    }

    SECTION("onivara_old_18, d1") {
        pmfe::ParameterVector params(mpq_class(-27, 10), mpq_class(-69, 10), 10, 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-3723, 10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(..........(...)..................................(...)...............).");
    }

    SECTION("onivara_old_19, d1") {
        pmfe::ParameterVector params(mpq_class(67, 10), 10, mpq_class(-21, 5), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-3873, 100));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((((((((.((....)))))(((....)))(.((....)))))((((.......)))))))))))).");
    }

    SECTION("onivara_old_20, d1") {
        pmfe::ParameterVector params(mpq_class(47, 10), 7, mpq_class(31, 5), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-108, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "..........(((.((....)))))(((((.......))))).((..(((((.......)))))..))....");
    }

    SECTION("onivara_old_21, d1") {
        pmfe::ParameterVector params(mpq_class(-73, 10), mpq_class(-3, 2), mpq_class(-9, 10), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-1939, 20));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(...(.........((....)).....(....)......)...(....).....................).");
    }

    SECTION("onivara_old_22, d1") {
        pmfe::ParameterVector params(mpq_class(29, 10), mpq_class(-1, 5), mpq_class(19, 10), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-1271, 50));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "((((((((......((....))...(((((.......))))).((((....))))........)))))))).");
    }

    SECTION("onivara_old_23, d1") {
        pmfe::ParameterVector params(mpq_class(42, 5), mpq_class(26, 5), mpq_class(15, 2), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-108, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "..........(((.((....)))))(((((.......))))).((..(((((.......)))))..))....");
    }

    SECTION("onivara_old_24, d1") {
        pmfe::ParameterVector params(mpq_class(63, 10), mpq_class(26, 5), mpq_class(17, 10), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-108, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "..........(((.((....)))))(((((.......))))).((..(((((.......)))))..))....");
    }

    SECTION("onivara_old_25, d1") {
        pmfe::ParameterVector params(9, mpq_class(-81, 10), mpq_class(-89, 10), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-12202, 25));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(..........(...)..................................(...)...............).");
    }

    SECTION("onivara_old_26, d1") {
        pmfe::ParameterVector params(mpq_class(1, 10), mpq_class(1, 5), mpq_class(16, 5), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-108, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "..........(((.((....)))))(((((.......))))).((..(((((.......)))))..))....");
    }

    SECTION("onivara_old_27, d1") {
        pmfe::ParameterVector params(-8, mpq_class(13, 5), mpq_class(38, 5), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-108, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "..........(((.((....)))))(((((.......))))).((..(((((.......)))))..))....");
    }

    SECTION("onivara_old_28, d1") {
        pmfe::ParameterVector params(mpq_class(-28, 5), mpq_class(11, 5), mpq_class(51, 10), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(108, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "..........(((.((....)))))(((((.......))))).((..(((((.......)))))..))....");
    }

    SECTION("onivara_old_29, d1") {
        pmfe::ParameterVector params(mpq_class(19, 10), 8, mpq_class(73, 10), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(108, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "..........(((.((....)))))(((((.......))))).((..(((((.......)))))..))....");
    }

    SECTION("onivara_old_30, d1") {
        pmfe::ParameterVector params(mpq_class(-43, 5), mpq_class(-91, 10), mpq_class(3, 2), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-5337/10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(..........(...)..................................(...)...............).");
    }

    SECTION("onivara_old_31, d1") {
        pmfe::ParameterVector params(mpq_class(3, 10), mpq_class(-13, 10), mpq_class(-9, 10), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-3901, 50));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(.............((....)).....................(....).....................).");
    }

    SECTION("onivara_old_32, d1") {
        pmfe::ParameterVector params(mpq_class(-7, 10), mpq_class(-1, 2), mpq_class(-3, 10), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-2297, 50));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "((((((((......((....))..((.(((.......)))...((((....))))....))..)))))))).");
    }

    SECTION("onivara_old_33, d1") {
        pmfe::ParameterVector params(0, -3, mpq_class(-2, 5), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-858, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(..............(....)......................(....).....................).");
    }

    SECTION("onivara_old_34, d1") {
        pmfe::ParameterVector params(mpq_class(6, 5), -1,  mpq_class(9, 10), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        mpq_class energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == mpq_class(-5629, 100));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "((((((((......((....)).....................(....)..............)))))))).");
    }
}
