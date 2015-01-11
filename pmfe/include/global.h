#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include "constants.h"
#include "pmfe_types.h"
#include <gmpxx.h>

namespace pmfe {
    using namespace std;

    extern unsigned char *RNA;
    extern int *structure;
    extern int* constraints;

    extern int g_nthreads;
    extern int g_unamode;
    extern int g_dangles;
    extern int g_mismatch;
    extern int g_verbose;
    extern int g_prefilter_mode;
    extern int g_prefilter1;
    extern int g_prefilter2;
    extern unsigned int chPairKey;

    extern int SHAPE_ENABLED;//0 means false and 1 means true
    extern int g_LIMIT_DISTANCE;
    extern int g_contactDistance;

    extern mpz_class INFINITY_;
    extern mpq_class inf;
    extern mpq_class zero;

    extern mpq_class dummy_scaling;
    extern mpq_class multiloop_penalty;
    extern mpq_class unpaired_penalty;
    extern mpq_class branch_penalty;

    // The possible base pairs are (A,U), (U,A), (C,G), (G,C), (G,U)
    //  and (U,G).
#define checkPair(i, j) (((((i)-(j)) % 2) == 1 || (((i)-(j)) % 2)== -1) && (!( ((i)==BASE_A && (j)==BASE_C) || ((i)==BASE_C && (j)==BASE_A) )))

    int canPair(int a, int b);
    void init_global_params(int len, ParameterVector params);
    void free_global_params();
    void print_sequence(int len);
    void print_structure(int len);
    void print_header() ;
    void print_gtfold_usage_help();

    int read_sequence_file(const char* filename, std::string& seq);
    bool encodeSequence(string seq);
    void save_ct_file(string outputFile, string seq, mpq_class energy) ;
    void save_ct_file(string outputFile, string seq, mpq_class energy, int *structure1);

    void init_checkPair();
    int  update_checkPair(int i, int j);
}

#endif
