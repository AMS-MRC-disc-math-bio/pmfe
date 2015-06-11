// Copyright (c) 2015 Andrew Gainer-Dewar

#include "catch.hpp"
#include <boost/filesystem.hpp>

#include "mfe.h"
#include "nndb_constants.h"
#include "nntm.h"
#include "pmfe_types.h"
#include "thread_pool.h"
#include "rational.h"

namespace fs = boost::filesystem;

// Build the thread pool
pmfe::SimpleThreadPool thread_pool;

TEST_CASE("A. tabira 5S MFE", "[mfe][biological][atabira][5S]") {
    // Load the sequence
    fs::path seqfile("test_seq/5S/a.tabira_5S.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 120);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-533, 10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((((.....((((((((((...((.((((......))))..))..)))...)))))))..(((((((...(.((..(..((....))..)..)))..))))))))))))))))..");
    }
}

TEST_CASE("C. diphtheriae tRNA MFE", "[mfe][biological][cdiphtheriae][tRNA]") {
    // Load the sequence
    fs::path seqfile("test_seq/tRNA/c.diphtheriae_tRNA.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 74);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-35));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "((((((((((.(.(.((((...((.(((.(((((((....)))).))).))).)))))).))))).))))))).");
    }
}

TEST_CASE("D. mobilis 5S MFE", "[mfe][biological][dmobilis][5S]") {
    // Load the sequence
    fs::path seqfile("test_seq/5S/d.mobilis_5S.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 133);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-463, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "..(((((((((((((((....((((((((..(((..((((..(((...)))..))))..)))...)))))))).((((((((((.((((((.((....)))))))).)))))))))))))))).)))))))))");
    }
}

TEST_CASE("E. coli 5S MFE", "[mfe][biological][ecoli][5S]") {
    // Load the sequence
    fs::path seqfile("test_seq/5S/e.coli_5S.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 120);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-549,10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "((((((((((.........((((....)))).((((((((((.....(((....)))...(((((.......))))).))))))))))..(((.(((....))))))..)))))))))).");
    }
}

TEST_CASE("G. arboreum 5S MFE", "[mfe][biological][garboreum][5S]") {
    // Load the sequence
    fs::path seqfile("test_seq/5S/g.arboreum_5S.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 120);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-407, 10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((((...((.((.(((....................(((((..(.((((((....)))))).).)))))...((((((.((....))))))))..))).)).)))))))))))..");
    }
}

TEST_CASE("H. sapiens tRNA MFE", "[mfe][biological][hsapiens][tRNA]") {
    // Load the sequence
    fs::path seqfile("test_seq/tRNA/h.sapiens_tRNA.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 72);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-263, 10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "((((.....)))).((((((.((.(((((((((..((((....))))..))).))))))))...))))))..");
    }
}

TEST_CASE("L. delbrueckii tRNA MFE", "[mfe][biological][ldelbrueckii][tRNA]") {
    // Load the sequence
    fs::path seqfile("test_seq/tRNA/l.delbrueckii_tRNA.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 72);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-239, 10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((...............(((.(((((.......))))).)))..(((.......)))..))))))).");
    }
}

TEST_CASE("O. nivara tRNA MFE", "[mfe][biological][onivara][tRNA]") {
    // Load the sequence
    fs::path seqfile("test_seq/tRNA/o.nivara_tRNA.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 73);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-271, 10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((..((((........))))((((((.......))))))....(((((.......)))))))))))).");
    }
}

TEST_CASE("R. norvegicus 5S MFE", "[mfe][biological][rnorvegicus][5S]") {
    // Load the sequence
    fs::path seqfile("test_seq/5S/r.norvegicus_5S.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 123);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-539, 10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((((....((((((((...((..((((..((....))..))))..))....)))))).))(((((((...(.((..(.(((....))).)..)))..)))))))))))))))).....");
    }
}

TEST_CASE("S. tokodaii tRNA MFE", "[mfe][biological][stokodaii][tRNA]") {
    // Load the sequence
    fs::path seqfile("test_seq/tRNA/s.tokodaii_tRNA.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 74);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-79, 2));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "((((((((.((((.((...)).))))..((((.(((.((((((.......)))))).))).)))))))))))).");
    }
}

TEST_CASE("A. suum 16S MFE", "[mfe][biological][asuum][16S]") {
    // Load the sequence
    fs::path seqfile("test_seq/16S/a.suum_16S.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 701);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-149));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "..((((.....))))...(((((((((((..............((((((..((.....))..))))))))))))))))).........((.(((((((....((((((..(((((..(((((.((((((((((((((......(((((.((((.((((......)))).)))).)))))(((..(((((..((((((...((.((((((....((((((((..(((((.........)))))..))))))))....)))))).)).......((((((........))))))..)))))).)))))..))).))))))))))).))).)))))..))))).....((((.((((((((((...(((..(((((..(((.((...((((....))))...))..(((((((.(((.((.((.((((....(((((((((.....((.(((..((((.((......)).))))..))).)))))))))))..)))).)))).)))))))))).....)))..))))......((((.((((.((.(((((((((....))))))).....)).)).)))).))))...................)..)))...))))))))))((..........))......)))).))))))....))))))).))..(.((((((((((....)))))))))).).....");
    }
}

TEST_CASE("B. bigemina 16S MFE", "[mfe][biological][bbigemina][16S]") {
    // Load the sequence
    fs::path seqfile("test_seq/16S/b.bigemina_16S.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 1701);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-573));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((......((((..(((((((((((((((.((.(((((.((.((((.....)))).)).))))).....((.((((((......(((((((((..((((((..(((..(((((((((..((((((.(((((......))))).))).))).....))))))).(((((((((......((..(((....)))..))..)))))))))................((((((((.(((.....))).))))))))..))..)))..)))))).)))))).))).....))))))))....)).))))))).(((((((((.(((((((((.....(((((((........((((..((((((((((........)))).)).).)))..))))(((.((.(((((....((((((...(((((...))))).))))))...................((((((..(((((.(((((((.((..((((((.(.....).))))))..)).))))))).))))))))))).))))).)).)))....))))))).(((.....))).........))))))))).(((......((((.(((........((((..(((......(((.......)))......))).))))........))).)))).......(((((((.(((((((....))).))))))))))).)))..((((((((((((((((((....))))))))))))).)))))..))))))))).((((((((..(((((......(((((((((.((((((.(((.((((......)))).)))..)))))).)))))....(((....)))..)))).....)))))..))))))))..(((....)))..)))))))).(((..((((((..(((((((.((((..(((..((((.....((((.((..(((((((..(((.......))).)))))))..))....((....))..))))....)))).))).)))).))))).))...))))))......(((((((..(.((((((((.((..(((((.....))))).)).((((((((.........((....))......(((((..((...))..)))))..))))))))..)))))))))..)))))))(((((((((((((.(((.((..((((((((..(.((.....((((.(((.....((....))))).))))....)))..))))))))..))..((((((..(((((((.(((..(((((.....((((..(((.....)))..))))..((((((((((...(((....)))..)).))))))))...)))))..))).((((((((....((((....)))).((.((..((...(((....)))..))..)).)).)))))))).(((((...(((((((.......)))))))....))))))))))))..........(((((((((.....)))))))))....))))))..)))...))))).))))))))....(((((.((....(((((((..(((..((....(...((..(((....((((....))))....)))..))...)....))..)))..)))))))...)).)))))..)))......))))(((((((((....)))))))))...)))...");
    }
}

TEST_CASE("C. elegans 16S MFE", "[mfe][biological][celegans][16S]") {
    // Load the sequence
    fs::path seqfile("test_seq/16S/c.elegans_16S.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 697);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-285, 2));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "...((((..(.....)..)))).......(((((.((((((.((((((..........)))))))))))))))))..........(((((....(((..(((.(((((((....((((....(((..(((((((((..(((((.(((..((..((....))..))..((((((.(((.((((((((((((((((((..........((((((((..........))))))))...))))))..)))))))))..))).))).))))))...))).)))))....))))...)))))..)))....)))).(((......(((((((....(((..((....))..)))...(((((.((.(((((..(((((..(((.((...((((....)))).......(((((.((...(((.((.(((((...((((((((((...(((((((((.........))))))))).....)))))))))).))))).)))))...)).)))))..)))))..)))..........(((((((.(((((((...(((((....)))))....(((((.....)))))..))))))).)))))))..))..))))).)).)))))..)))))))....)))...)))))))...)))..)))..))))).....(.((((((((((....)))))))))).)....");
    }
}

TEST_CASE("E. cuniculu 16S MFE", "[mfe][biological][ecuniculi][16S]") {
    // Load the sequence
    fs::path seqfile("test_seq/16S/e.cuniculi_16S.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 1295);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-487));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "........(((((((((.(((.((...(((.((((((((((((((((.....(((((((.((((((..(..((((....((((((((.((((((....)))..))).)))))))).....((.((.((.((((((....)))))).)).)).))..))))..))))))).)))))))..(((((..((..((.((...((((((((.......))))))))...)).))..)).)))))........))))))))))))))))((.((((((...(((((.(((((((((((.(.(..((.(((.(((((..((((((((..(((....))))))).........))))..))))).)))..(..((((((.((((((...(((((..(((((..(..(.(((.((((.((((.(((((((...........)))).))).)))).))))))).)..)..))))).)))))...)))).)).))))))..)............))..).).)))))))))))..(((((((.((((((((.((((((((.(((...((.(.(....).).))..(((....((((.((.((.((.(((((((((..(((..((((((.(((((((.(((.(......))))...((.((...(((((....)))))...))))...))))))).)).))))......)))...(((((((.((.....)).)))))))......((((...((((((...(((((..((((((..((((....))))...))).......)))..))))).....((((.....)))).(((((((((((.......((((.(.(((((((.....((.((((((.((((((((...((((.((((.(((..(....(((((.(((...))).)))))......)..))).)))).)))).((((((.....(....)...))))))........)))))))).))))))..)).(((..((........))..)))....))))))).).))))..)))))))))))))))))...))))..))))))))).)))))).))))..((((((.(....).))))))....))).....))).......(((((((.....)))))))......)))))))))))))))))))))))...)))))..(((((.(..((((((.(.((((..(((((....)))))..)))).).))))))..).)))))..........)).)))).))..)))....)).))))))))))))...");
    }
}

TEST_CASE("E. hexamita 16S MFE", "[mfe][biological][ehexamita][16S]") {
    // Load the sequence
    fs::path seqfile("test_seq/16S/e.hexamita_16S.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 1550);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-2824, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "...(((((.......)))))...(((..(((((((((........((((.((((((((..((((.((((..(((....)))..(((((((((((..((((...((....))))))..)))).......))))))).)))).))))(((((((((.((.....((((((......))).(((((((...((((((...(((((.(.....(((.(((.(((.(((.....(((((((((((...((((.(((((((.((....(((((((((.....(((((((....(((((..(((...((..((((((.(((......))).....((((((.((.....))))))))))))))..))..(((.((.(((((.(..((((((...(((((....))))).))))))..)............((((((....)))))).....(((..((((((.(......)))))))...(((.....)))....)))...))))).)).))))))..))))).......))))))))))))))))..)).))).))))..)))))))).)))...(((((((((((((..((((((((.....)))))))).(((((...((....))((((.......))))....)))))...(((((((((....))))))))).((.......((((((((...((((((.((.(((((((((...................))))))))).))........((....)).))))))..))))))))..(((((((((((((....((((((((..(((((............)))))..))))))))......(((((((.....(((.(((((((....))))))).)))...((....)).((((((((((.(..((.((.........)).))..)))))))))))..)))))))......((((.((((.....)))).))))...))))...))))).))))............)).(.((((((((.((((.((.((((((((((....(((....((((((.((((((.((((..(((((.((.((.....))..)).))))).)))).)))))).))))))..)))..((.((((((((((((((.....)))))...(((((((.....))))).))....)))))).))).)).)))))))..))))).)))))).)))))).)..((((...((((.((((((((..((((.((..(((.(((....)))..)))....))))))....)))))))).)))).((((.......(((...(((((((....)))))))..)))....)))).)))).....)))))))))))))..))))...))).))).))))))).))))))))))))))))))....))).....)).)))))))))..............)))))))).((..((((.(...(((((.....))).))..).)))).)).......))))((((((((....))))))))..)))).))))))))");
    }
}

TEST_CASE("G. ardaea 16S MFE", "[mfe][biological][gardaea][16S]") {
    // Load the sequence
    fs::path seqfile("test_seq/16S/g.ardaea_16S.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 1435);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-7619, 10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "..(((..(.((.((((((.((((((..((..((((((..(....)..(((.(((.(((....)))))).)))...))))))..))..........((..((.(.(((.(((.(((((((..(((..(((((((..((((((....))))))....(((..((((....)))).((.(.((.((((.((((((((........((((.((((..(((((((((.....(((.((((((((.......)))))))))))......((.(((((..(((((....((((((((((((..(((((((........)))))))...)))......((..(((((.......))))))).(((((.(((((.....((.((.(((...((((((((...))))))))..))).))))...))))).)))))))))))))).....)))))..))))))).))))))).))..))))..))))..)))))))))))).)).).)).))).....)))))))))).))))))).))).))).).))..))..))))))..)))))).)))..)))..(((.((.((((.(((.(.(..(((.((((((((((..((((((((...(((((((((.................))))....)))))(((((.((((((((((..(((((..((((((((..(((.(((((..((((((((..(.((((((.((..(.((((((((....((((.(((....(((((((.((((((.(((((.((((..(((..((....))...)))..).))).))))).)))))).((((.(((((....)))))..))))...)))))))(((((((.(((((..(((((.....)))))(((((....(((((((..((((..(.(((....)))...)..))))..)))..)))).)))))))))).)))))))...))).))))))))))))))).))).))).).((.((((((....)))(((.((......))))).))).))....))))))))...((((....))))..)))).).)))..))))))))...(((((.(.......).)))))...))))))))))).)))).))))).)).))))))(((.(((((((...(((.((((.(((((....))))).))).).)))....)))))))))).(((((.(((((((..((....))))))))).)))))..(((((....(((((((((.....)))))))))))))).)))))))))))))...).).)))))))...(((((.((..(..(((((((.(((..(((.(((...((.....))...))).)))..))).)))))))..).)).)))))..)).)))...(((((((((((((....))))))))))....)))..");
    }
}

TEST_CASE("G. intestinalis 16S MFE", "[mfe][biological][gintestinalis][16S]") {
    // Load the sequence
    fs::path seqfile("test_seq/16S/g.intestinalis_16S.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 1452);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-8027, 10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                ".....((..(((((.(((...((((.(..(((((((((..(.((.((((......((((((((...((((..(((((((.(((((((((...((.(((((.......))))).)).))))..)))))....))))).))..))))((((((.....)).))))..)))))))).(((..(((((.((.(((.(((((.(((.(((((.(((((..((((((((((.....((....))....(((((((((((((..(((((....((((((.(..((((....))))..(((((((..(((((((........)))))))...)))..((.((((....((((((((.(((((((..(((((...))))).(((..((.(((....(((....))).((.((((.....)))))).)))))..)))..))))))).))).)))))..((((..(((......(((((((((.((((((....)))))))))).)))))......)))..)))).......))))))....(....).))))..).)))))).)))))...((.(((((.(((((....))))).)))))))))))))))))))).)))))).))))..)))))((((((((((....((((((((((..((((.....(..(((((((....((.((.((((((..((((....))))......((((..((.....))..))))..)))))).))))....)).)))))..)..))))...))))))))))....))))))).)))...))))).)))..))))).)))..))))))).((....))))).)))).))..)..)))))))))..).).)))..))).....((((..(((...((((((((((((....((((....((((.....))))..)))).......(((((..(((....))).((((.(((((((.((.......((((((..((((((((((((.(.(((((((((.((((((((..((.....)))))).))))..))))))))).).)))))))....)))))....)))))).(((((((((((((((((.(((((....)))))...(((((((((........(.((((((.((.((((....(((...)))....((((((((...(((((...)))))..)))))))).)))).)))))))).))))))))))..((.((..((((((..((....))))))))..)).))..)))))))).)))))))))..)).))))))))))).)))))...)))))))............))))).(((((.((..(..(((((((.(((..((((((((......))..))).)))..))).)))))))..).)).)))))..)))....))))....((((((((....)))))))).)))))))..");
    }
}

TEST_CASE("G. muris 16S MFE", "[mfe][biological][gmuris][16S]") {
    // Load the sequence
    fs::path seqfile("test_seq/16S/g.muris_16S.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 1432);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-582));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((....((((((..(((((((..((((((((..((......(((((((...((....))..)).....)))))....((((((.(((....)))))))))...)).)))))...))).....((....((((((....))))))))......)))))))....)))))).......((((((......)))..))))))))))....(((((((((((((((..((.(((((.......))))).)).))).))))))))..))))(((((...((((((((((.....(((((((........)))))))...((((((...((((.(((((.((....(((((..((((......))))..)))))...)).....((((.(((..((......((((((.........))))))..))..))).))))...))))).))))......((((((......((((((((((.((((....)))))))))).))))....(((((....(((...(((((..(((((((.....(((((((.....(((((.((((((((((....((((((((.((..(((........)))..(((((((((.....(((((.(((...................))).))))).((((((...(..((..((.((..((....((((..........(((((...(((((.(.((((..((((((.((..(.((((((((.......((((.....(((((((.((.((((.((((((.((((.....))))...((....))......))).))).)))).))..((((.(((((....)))))..))))...)))))))..((..((((.......(.(((.(((((.............))))).))))....))))..))..))))...(((((((.....(((.((((.....)))).)))...)))))))..))))))))))).))).)))...)))).)))))).....(((((((((..........))))).)))).....)))))))))..((.(((((.....))))).))..))..)))).))..)..)))))).)))))))))..)).)))))))).....)))))).))))((((((.(((...(.(((((((((....))))))))).).)))..))))))..(((.((((.((((.((((((...(((..(....)..)))...)))))).(((((((((.....)))))))))...)))))))).)))...........)))))...)))))))..((.(..((...))..)))..)))))))...)))))...)))....))))).))))))...))))))....)))...)))).)))....((((((((((....)))))))))))))))....");
    }
}

TEST_CASE("H. volcanii 16S MFE", "[mfe][biological][hvolcanii][16S]") {
    // Load the sequence
    fs::path seqfile("test_seq/16S/h.volcanii_16S.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 1474);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-687));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                ".....(((.(((((.......((((((..((.(((((((((......(((.((((..((((((((((....)))))).))))....(((.......((..((((......(((((((.((((..((.((((((....)))))).))...))))....(((((((..((......)).))))))).....(((....))))))))))...))))..))(((((((((((..((((((((.......)))))))))))..)))))))).((((((((....)))...)))))))).((((((((......)))))))).((((....))))..)))))))......((((.....(((((....)))))..)))))))))))))..((((((((.....)))))))).(((((((((((.....)))))))))))...((((((.......((((((...((((...))))..)))))))))))).))..))))))...........((((....(((((.(.((((((((((..(((((((((((...((((((((.....))))))))..))))))))..))).))))))...((((((.....((((((((((.((((((((((.(((......)))....)))))))))).))).......((....)).)))))))....))).)))..)))).).)))))..((((((((.((....((((.........))))..))))))))))....((((((((((((.((((..((((((((.....))))))))..))))...((....))....)))))))..(((...((((.(((((....)))))..))))..)))....((((((((((((....(((((.((....)).))))).....((((....(((((.((.(((..(((((((...((((((((.((..(((.(((((...))))).))).)).))))))))...((.((((..((((((((((..((((((.((....)).(.((((......))))).))))))..(.((((..(((((((((((((((..((((....))))))).)))).))))))))..((.(((((.....))))).))......))))).))))).)))))...))))))....))))))).......(((((((((.(((((((((.....(((((((..((((....))))..)).)))))......))))))))).((((((..(((((((........))))))).....))))))...)))))..))))........))))).))))).....)))).)))))))........))))).....((((.(((..(((((((((((((..(((((((....)))))))..)))))))))))))..))).))))..)))))....)))).((((((((....)))))))).)))))))).....");
    }
}

TEST_CASE("V. necatrix 16S MFE", "[mfe][biological][vnecatrix][16S]") {
    // Load the sequence
    fs::path seqfile("test_seq/16S/v.necatrix_16S.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 1244);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-1628, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "........((((((((((((((((...((((....(((.........(((.(((((((((((.((((((......(((((((((((.((((.((((((.(((((....)))))..(((((.((((((....(((....(((....(((((((..((...((((.(((..((((((..((((((((.......)))))))))))..)))..)))...(((((((((((.((....(((......)))((((.(((.....(((.((.(((.(((..((((.......))))..)))))).)))))......))).)))).....)).)))))))))))........))))....))...)))))))..)))....))).(..((((.......))))..)...)))))).))))).......)))))).))))....))))).))))))......)))))).......(((((((((..((((.((((((((.((.((((..((.....................))..)))).)).......(((....)))))))))))))))...))))))))))))))))))..))))))))))))...))).)))).......((...((((((((((((.(((((...(((((.(((((...)))))..((((((((((...........(((((....)))))...(((((((..((((....))))................(((((.......((....)).((((...(((((((((((..((.((....((..(((....(((((.((((((((((..(((.(.....)..)))..)))))).)))).)))))....)))......(((((((.((((..............((((((......)))))).((((((.((..((.((.((.....)).)).))..)).))))))......)))))))))))..))....)).))..)))))))))))......(((((..........((((((((..(.......)..)).))))))....((((..(((..(((.....(((((((((.....)))))))))....))))))..)))).)))))....))))))))))))))))....(((((((((((((..((.((((....)))).))..))))))))))))).)))))))))).)))))..)))))))))))))))))...))...)))))))))...");
    }
}

TEST_CASE("Z. mays 16S MFE", "[mfe][biological][zmays][16S]") {
    // Load the sequence
    fs::path seqfile("test_seq/16S/z.mays_16S.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 1962);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-7419, 10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "....((((((((.......)))))))).((((((((.((.....))....(((..(((((.(((..(((((((((((((.(((((....((((.......(((((..(((...((((((.(((((((....))))))).).)))))..............((((((((((((((((...((((.((((((...((((.((((((..(.(.(((((..(.((((((..((....))..((((((((((..((((.(((((.((((.....((((.((((...))))..)))).(((((((((..((((((((.......)))))))))))..))))))..))))..)))))..))))..((((((......)).)))).....)))))))))))))))).)..)).))).).).......((((...(((((....))))).))))((....((((....))))......))..))))))))))...(((..(((.((((((.(............).)))))).)))..)))..))))))..))))..)))))).(((.(((...((((.(.((((...((((((((.......((.(((((..(((.(((......(((((((.((.((((.((((....((((((.....))))))..)))).))))..)).)))))))...))).(((...(((((((((..((((((((.((.(((.....))).)).))..)))))).....))....((....)).)))))))..))))))....))))).)).(((....))).))))))))..((((.........))))...)))).).)))).......(((((((.(((((((((((((.((.((...)))).)))))))))))))...((....))....))))))).)))))).((((.(((((....)))))..)))).........))))))))))......)))...))))).(((((((...........))).))))...)))).))))))))))).....)))))))....(((((.....)))))......(((...((((((((((((((....(((.(((.(((((.((((((....)))....((......))....))).(((((((....)))))))...)))))))).....(((((((.(((((...((..((((..(((.((((((..(((((.(.........))))))..)))))).)))..))))..)).....))))).)).))))).)))....)))))))....)))))))...)))))).)))))...........((((((((((((((((..(((((........(((((((.((((((((.(((((((((..((((.(...((((...))))....((((((.(((.....))).))....(((((((...(((((.((((((..(((..(((((....((((..(((..(((((...............(((.(((((.((((((..((((.((.((((((...)))))).))...(((.((((........)))).)))..))))..)))).)))).))).))).((((.....)))).........)).)))..)))......((((((..((((((((..((.....))..))))))))...))))))..((((((.(((((...(((......)))...(((((((.........)))))))...))))).))))))....))))..))))).))))))))))))))..)))))))))))).))))....)))))..))))..))))).)))))).)))))))))((((((..((((((......))))))..)))))).......))).)))))))))..))))))).)).)))..(((..((((((((((....))))))))))..)))..)))...");
    }
}

TEST_CASE("Combinatorial sequence MFE", "[mfe][synthetic][combinatorial]") {
    // Load the sequence
    fs::path seqfile("test_seq/synthetic/test_combinatorial.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 60);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-149, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "((((....((((....))))....((((....((((....))))....))))....))))");
    }
}

TEST_CASE("Randomly generated sequence MFE", "[mfe][synthetic][random]") {
    // Load the sequence
    fs::path seqfile("test_seq/synthetic/test_random.fasta");
    pmfe::RNASequence seq(seqfile);

    // Some basic sanity checks
    REQUIRE(seq.len() == 75);

    SECTION("Turner99 published parameters") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-81, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                ".(((.(((....(((....((.....))..))).....))))))................((......)).....");
    }
}

TEST_CASE("O. nivara tRNA (old) MFE", "[mfe][biological][onivara-old][trna]") {
    // Load the sequence
    std::string theseq = "AUCAGAGUGGCGCAGCGGAAGCGUGGUGGGCCCAUAACCCACAGGUCCCAGGAUCGAAACCUGGCUCUGAUA";
    pmfe::RNASequence seq(theseq);

    // Some basic sanity checks
    REQUIRE(seq.len() == 72);

    SECTION("onivara_old_classical, d1") {
        pmfe::Turner99 constants(thread_pool);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-283, 10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "((((((((..(((.((....)))))(((((.......))))).((((....))))........)))))))).");
    }

    SECTION("onivara_old_1, d1") {
        pmfe::ParameterVector params(1, 1, 1, 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-239, 10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((.(((.((((....)).))(((((.......)))))..)))(((((.......)))))))))))).");
    }

    SECTION("onivara_old_2, d1") {
        pmfe::ParameterVector params(-1, 1, 1, 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-697, 25));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((.(((.((((....)).))(((((.......)))))..)))(((((.......)))))))))))).");
    }

    SECTION("onivara_old_3, d1") {
        pmfe::ParameterVector params(-1, 0, 1, 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-1669, 50));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((.(((...((....))...(((((.......)))))..)))(((((.......)))))))))))).");
    }

    SECTION("onivara_old_4, d1") {
        pmfe::ParameterVector params(1, 0, -1, 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-2067, 50));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((.(((...((....))...(((((.......)))))..)))(((((.......)))))))))))).");
    }

    SECTION("onivara_old_5, d1") {
        pmfe::ParameterVector params(pmfe::Rational(1, 10), 1, pmfe::Rational(-1, 10), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-806, 25));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((.(((.((((....)).))(((((.......)))))..)))(((((.......)))))))))))).");
    }

    SECTION("onivara_old_6, d1") {
        pmfe::ParameterVector params(3, pmfe::Rational(-1, 10), 1, 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-2729, 100));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "((((((((..(((.((....)))))(((((.......))))).((((....))))........)))))))).");
    }

    SECTION("onivara_old_7, d1") {
        pmfe::ParameterVector params(2, pmfe::Rational(1, 10), 5, 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-108, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "..........(((.((....)))))(((((.......))))).((..(((((.......)))))..))....");
    }

    SECTION("onivara_old_8, d1") {
        pmfe::ParameterVector params(2, pmfe::Rational(1, 10), 0, 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-65, 2));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((.(((...((....))...(((((.......)))))..)))(((((.......)))))))))))).");
    }

    SECTION("onivara_old_9, d1") {
        pmfe::ParameterVector params(0, pmfe::Rational(1, 10), 0, 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-73, 2));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((.(((...((....))...(((((.......)))))..)))(((((.......)))))))))))).");
    }

    SECTION("onivara_old_10, d1") {
        pmfe::ParameterVector params(pmfe::Rational(-1, 10), pmfe::Rational(1, 10), 2, 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-617, 25));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((.(((...((....))...(((((.......)))))..)))(((((.......)))))))))))).");
    }

    SECTION("onivara_old_11, d1") {
        pmfe::ParameterVector params(pmfe::Rational(17, 5), pmfe::Rational(1, 10), pmfe::Rational(2, 5), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-55, 2));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((((((((.((....)))))(((((.......))))).((((....))))....))..)))))))).");
    }


    SECTION("onivara_old_12, d1") {
        pmfe::ParameterVector params(pmfe::Rational(17, 5), pmfe::Rational(1, 5), pmfe::Rational(2, 5), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-27, 1));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((((((((.((....)))))(((((.......))))).((((....))))....))..)))))))).");
    }

    SECTION("onivara_old_13, d1") {
        pmfe::ParameterVector params(pmfe::Rational(16, 5), pmfe::Rational(1, 5), pmfe::Rational(2, 5), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-136, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((((((((.((....)))))(((((.......))))).((((....))))....))..)))))))).");
    }

    SECTION("onivara_old_14, d1") {
        pmfe::ParameterVector params(pmfe::Rational(16, 5), 0, pmfe::Rational(2, 5), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-5034, 25));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "((....)((((...)(....))((...)((.......)((...)))))(((...)(....))(....)))).");
    }

    SECTION("onivara_old_15, d1") {
        pmfe::ParameterVector params(pmfe::Rational(43, 5), pmfe::Rational(31, 10), pmfe::Rational(-131, 10), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-143, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((.(((...((....))...(((((.......)))))..)))(((((.......)))))))))))).");
    }

    SECTION("onivara_old_16, d1") {
        pmfe::ParameterVector params(pmfe::Rational(18, 5), pmfe::Rational(-63, 10), 2, 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-35461, 100));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(..........(...)..................................(...)...............).");
    }

    SECTION("onivara_17, d1") {
        pmfe::ParameterVector params(pmfe::Rational(-7, 5), pmfe::Rational(39, 10), -10, 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-5172, 25));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "((....)((((...)(....))((...)((.......)((...)))))(((...)(....))(....)))).");
    }

    SECTION("onivara_old_18, d1") {
        pmfe::ParameterVector params(pmfe::Rational(-27, 10), pmfe::Rational(-69, 10), 10, 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-3723, 10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(..........(...)..................................(...)...............).");
    }

    SECTION("onivara_old_19, d1") {
        pmfe::ParameterVector params(pmfe::Rational(67, 10), 10, pmfe::Rational(-21, 5), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-3873, 100));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(((((((((((((.((....)))))(((....)))(.((....)))))((((.......)))))))))))).");
    }

    SECTION("onivara_old_20, d1") {
        pmfe::ParameterVector params(pmfe::Rational(47, 10), 7, pmfe::Rational(31, 5), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-108, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "..........(((.((....)))))(((((.......))))).((..(((((.......)))))..))....");
    }

    SECTION("onivara_old_21, d1") {
        pmfe::ParameterVector params(pmfe::Rational(-73, 10), pmfe::Rational(-3, 2), pmfe::Rational(-9, 10), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-1939, 20));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(...(.........((....)).....(....)......)...(....).....................).");
    }

    SECTION("onivara_old_22, d1") {
        pmfe::ParameterVector params(pmfe::Rational(29, 10), pmfe::Rational(-1, 5), pmfe::Rational(19, 10), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-1271, 50));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "((((((((......((....))...(((((.......))))).((((....))))........)))))))).");
    }

    SECTION("onivara_old_23, d1") {
        pmfe::ParameterVector params(pmfe::Rational(42, 5), pmfe::Rational(26, 5), pmfe::Rational(15, 2), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-108, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "..........(((.((....)))))(((((.......))))).((..(((((.......)))))..))....");
    }

    SECTION("onivara_old_24, d1") {
        pmfe::ParameterVector params(pmfe::Rational(63, 10), pmfe::Rational(26, 5), pmfe::Rational(17, 10), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-108, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "..........(((.((....)))))(((((.......))))).((..(((((.......)))))..))....");
    }

    SECTION("onivara_old_25, d1") {
        pmfe::ParameterVector params(9, pmfe::Rational(-81, 10), pmfe::Rational(-89, 10), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-12202, 25));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(..........(...)..................................(...)...............).");
    }

    SECTION("onivara_old_26, d1") {
        pmfe::ParameterVector params(pmfe::Rational(1, 10), pmfe::Rational(1, 5), pmfe::Rational(16, 5), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-108, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "..........(((.((....)))))(((((.......))))).((..(((((.......)))))..))....");
    }

    SECTION("onivara_old_27, d1") {
        pmfe::ParameterVector params(-8, pmfe::Rational(13, 5), pmfe::Rational(38, 5), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-108, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "..........(((.((....)))))(((((.......))))).((..(((((.......)))))..))....");
    }

    SECTION("onivara_old_28, d1") {
        pmfe::ParameterVector params(pmfe::Rational(-28, 5), pmfe::Rational(11, 5), pmfe::Rational(51, 10), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(108, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "..........(((.((....)))))(((((.......))))).((..(((((.......)))))..))....");
    }

    SECTION("onivara_old_29, d1") {
        pmfe::ParameterVector params(pmfe::Rational(19, 10), 8, pmfe::Rational(73, 10), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(108, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "..........(((.((....)))))(((((.......))))).((..(((((.......)))))..))....");
    }

    SECTION("onivara_old_30, d1") {
        pmfe::ParameterVector params(pmfe::Rational(-43, 5), pmfe::Rational(-91, 10), pmfe::Rational(3, 2), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-5337/10));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(..........(...)..................................(...)...............).");
    }

    SECTION("onivara_old_31, d1") {
        pmfe::ParameterVector params(pmfe::Rational(3, 10), pmfe::Rational(-13, 10), pmfe::Rational(-9, 10), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-3901, 50));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(.............((....)).....................(....).....................).");
    }

    SECTION("onivara_old_32, d1") {
        pmfe::ParameterVector params(pmfe::Rational(-7, 10), pmfe::Rational(-1, 2), pmfe::Rational(-3, 10), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-2297, 50));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "((((((((......((....))..((.(((.......)))...((((....))))....))..)))))))).");
    }

    SECTION("onivara_old_33, d1") {
        pmfe::ParameterVector params(0, -3, pmfe::Rational(-2, 5), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-858, 5));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "(..............(....)......................(....).....................).");
    }

    SECTION("onivara_old_34, d1") {
        pmfe::ParameterVector params(pmfe::Rational(6, 5), -1,  pmfe::Rational(9, 10), 1);
        pmfe::Turner99 constants(thread_pool, params);
        pmfe::NNTM energy_model(constants, pmfe::CHOOSE_DANGLE, thread_pool);

        pmfe::RNASequenceWithTables seq_annotated = energy_model.energy_tables(seq);

        pmfe::Rational energy = energy_model.minimum_energy(seq_annotated);

        REQUIRE(energy == pmfe::Rational(-5629, 100));

        pmfe::RNAStructureWithScore scored_structure = energy_model.mfe_structure(seq_annotated);

        REQUIRE(scored_structure.old_string() ==
                "((((((((......((....)).....................(....)..............)))))))).");
    }
}
