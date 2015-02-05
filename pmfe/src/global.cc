#include <stdio.h>
#include <stdlib.h>
#include <fstream>

#include "utils.h"
#include "global.h"
#include "constraints.h"
#include "pmfe_types.h"
#include <gmpxx.h>

namespace pmfe {
    unsigned char *RNA;
    int *structure;
    unsigned int chPairKey;

    int g_nthreads;
    dangle_mode g_dangles;
    int g_unamode;
    int g_mismatch;
    int g_verbose = 0;
    int g_prefilter_mode;
    int g_prefilter1;
    int g_prefilter2;

    int SHAPE_ENABLED = 0;
    int g_LIMIT_DISTANCE;
    int g_contactDistance;
    int g_bignumprecision = 512;

    mpz_class INFINITY_;

    void init_global_params(int len, ParameterVector params) {
        RNA = new unsigned char[len+1];
        if (RNA == NULL) {
            perror("Cannot allocate variable 'RNA'");
            exit(-1);
        }
        structure = new int[len+1];
        if (structure == NULL) {
            perror("Cannot allocate variable 'structure'");
            exit(-1);
        }

        dummy_scaling = params.dummy_scaling;
        multiloop_penalty = params.multiloop_penalty;
        unpaired_penalty = params.unpaired_penalty;
        branch_penalty = params.branch_penalty;

        INFINITY_ = mpz_class("9999999999999");

        if (abs(params.dummy_scaling > 1)) {
            INFINITY_ *= abs(params.dummy_scaling);
        }

        init_checkPair();
    }

    void free_global_params() {
        delete[] structure;
        delete[] RNA;
    }

    void print_sequence(int len) {
        int i;
        for (i = 1; i <= len; i++) {
            if (RNA[i] == BASE_A)
                printf("A");
            else if (RNA[i] == BASE_C)
                printf("C");
            else if (RNA[i] == BASE_G)
                printf("G");
            else if (RNA[i] == BASE_U)
                printf("U");
            else
                printf("N");
        }
        printf("\n");
    }

    void print_structure(int len) {
        int i = 1;
        for (i = 1; i <= len; i++)
        {
            if (structure[i] > 0 && structure[i] > i)
                printf("(");
            else if (structure[i] > 0 && structure[i] < i)
                printf(")");
            else
                printf(".");
        }
        printf("\n");
    }

    int read_sequence_file(const char* filename, std::string& seq) {
        seq = "";

        ifstream fs;
        fs.open(filename, ios::in);
        if (!fs.good()) return FAILURE;

        string line;
        while(fs.good()) {
            getline(fs, line);
            // exclude lines starting with FASTA comment characters
            if(line[0] != ';' && line[0] != '>' && line.length() > 0)
                seq += line;
        }

        fs.close();

        size_t loc;
        while((loc = seq.find(" ")) != string::npos)
            seq.erase(loc, 1);

        return SUCCESS;
    }

    bool encodeSequence(string seq) {
        unsigned int unspecifiedBaseCount=0;
        unsigned int unspecifiedBases[seq.length()];

        for(unsigned int i=1; i<=seq.length(); i++) {
            RNA[i] = encode(seq[i-1]);

            // die on non-IUPAC codes
            if (RNA[i]=='X') {
                fprintf(stderr, "ERROR: Non-IUPAC nucleotide code at position: %d (%c)\n", i, seq[i-1]);
                fprintf(stderr, "See http://www.bioinformatics.org/sms/iupac.html for valid codes.\n");
                return false;
            }

            // add non-canonical IUPAC codes to the warning list
            if(!isWatsonCrickBase(seq[i-1]))
                unspecifiedBases[unspecifiedBaseCount++]=i;
        }

        // just print a warning for non-canonical IUPAC codes
        if(unspecifiedBaseCount > 0) {
            printf("\nIncompletely-specified IUPAC codes have been detected at position%s: ", unspecifiedBaseCount == 1 ? "" : "s");

            // put [0] first for nice comma separation
            printf("%d (%c)", unspecifiedBases[0], seq.at(unspecifiedBases[0]-1));
            for(unsigned int i=1; i<unspecifiedBaseCount; i++)
                printf(", %d (%c)", unspecifiedBases[i], seq[unspecifiedBases[i]-1]);

            printf("\nPlease replace with fully-specified IUPAC codes (A,C,G,U,T) and retry.\n");
            return false;
        }

        return true;
    }


    void print_header() {
        printf("GTfold: A Scalable Multicore Code for RNA Secondary Structure Prediction\n");
        printf("(c) 2007-2011  D.A. Bader, C.E. Heitsch, S.C. Harvey\n");
        printf("Georgia Institute of Technology\n\n");
    }

    void print_gtfold_usage_help() {
        //printf("There are three programs provided by GTfold. \n 1. gtmfe: Program implementing and supporting various options for calculating minimum free energy and MFE structure for a given sequence. Run 'gtmfe --help' for more help on this.\n 2. gtboltzmann: Program implementing and supporting various options and features for sampling structures for a sequence stochastically with help of partition function. Run 'gtboltzmann --help' for more help on this.\n 3. gtsubopt: Program implementing and supporting vcarious options and features for getting sub-optimal structures for a given sequence. Run 'gtsubopt --help' for more help on this.\n\nPlease visit our sourceforge page for more details,\nhttp://gtfold.sourceforge.net/\n\n");
        printf("The GTfold package consists of three programs.\n 1. gtmfe: Computes MFE structure for a given sequence. Use the command ‘gtmfe --help’ for more details.\n 2. gtboltzmann: Supports various features such as stochastic sampling and computation of base pair probabilities under the Boltzmann Distribution. Use the command ‘gtboltzmann --help’ for more details\n 3. gtsubopt: Computes all structures with in a given energy range of the MFE for a given sequence. Run ‘gtsubopt --help’ for more details.\n\n\n Please visit our sourceforge page for full details on all three programs\n http://gtfold.sourceforge.net/\n");
    }


    void save_ct_file(string outputFile, string seq, mpq_class energy) {

        ofstream outfile;
        outfile.open(outputFile.c_str());

        outfile << seq.length() << "\t  dG = " << energy.get_str() << endl;
        //outfile << seq.length() << "\tdG = " << energy/100.0 << "\t" << seqfile << endl;

        unsigned int i = 1;
        for(i=1; i <= seq.length(); i++)
            outfile << i << "\t" << seq[i-1] << "\t" << i-1 << "\t" << (i+1)%(seq.length()+1) << "\t" << structure[i] << "\t" << i << endl;

        outfile.close();
    }

    void save_ct_file(string outputFile, string seq, mpq_class energy, int *structure1) {
        ofstream outfile;
        outfile.open(outputFile.c_str());

        outfile << seq.length() << "\t  dG = " << energy.get_str() << endl;
        //outfile << seq.length() << "\tdG = " << energy/100.0 << "\t" << seqfile << endl;

        unsigned int i = 1;
        for(i=1; i <= seq.length(); i++)
            outfile << i << "\t" << seq[i-1] << "\t" << i-1 << "\t" << (i+1)%(seq.length()+1) << "\t" << structure1[i] << "\t" << i << endl;

        outfile.close();
    }

    void init_checkPair() {
        int i, j;
        chPairKey = 0;
        for (i = 0; i < 4; i++)
            for (j = 0; j < 4; j++)
                chPairKey += checkPair(i, j) << ((i << 2) + j);
    }

    int update_checkPair(int i, int j) {
        int r = 0;
        if (!((i >= 0 && i <=3 )&&(j >=0 && j <=3)))
            return r;
        if (!(chPairKey & (1 << ((i << 2) + j)))) {
            chPairKey += 1 << ((i << 2) + j);
            r = 1;
        }
        return r;
    }

    int canPair(int a, int b) {
        return (chPairKey & (1 << (((a)<<2) + (b))));
    }

    dangle_mode convert_to_dangle_mode(int n) {
        switch (n) {
        case 0:
            return NO_DANGLE;
            break;

        case 1:
            return CHOOSE_DANGLE;
            break;

        case 2:
            return BOTH_DANGLE;
            break;

        default:
            exit(EXIT_FAILURE);
            break;
        }
    }
}
