#include "StructureReader.h"
#include "RNAScoring.h"
#include "TreeScoring.h"
//#include <ctype.h>
//#include <limits.h>
#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
#include "loader.h"
//#include "rna-data.h"
//#include "constants.h"
#include <cstdlib>
#include <iostream>
#include<string.h>
#include "LoopScoring.h"
#include "options.h"
#include <math.h>

#include <string>
#include <vector>
#include <gmpxx.h>

using namespace std;
using namespace rnascoring;

namespace rnascoring
{
    extern char paramDir[200];

    mpq_class get_classical_score(std::string structfile, std::string paramdir, int dangle_model) {
        // Force parameters for iB4e use

        switch (dangle_model) {
        case 0:
            NODANGLEMODE=1;
            break;
        case 2:
            D2MODE=1;
            break;
        default: // includes dangle_model=1
            DEFAULTMODE=1;
            break;
        }

        // This code uses char arrays, and I don't feel like rewriting the whole thing
        // to use C++ strings
        std::vector<char> structfile_char (structfile.begin(), structfile.end());
        structfile_char.push_back('\0');
        std::vector<char> paramdir_char (paramdir.begin(), paramdir.end());
        paramdir_char.push_back('\0');

        ResultBundle* resultBundle = CreateFromFile(&*structfile_char.begin());
        int length = resultBundle->length;
        TreeNode* tree = resultBundle->treenode;
        int* RNA = resultBundle->RNA_seq;

        nndb_constants* param = populate(&*paramdir_char.begin(), 1);

        int tree_score_fixed = ScoreNode(tree, resultBundle->RNA_seq, param, length);
        mpq_class tree_score = mpq_class(tree_score_fixed, 100);
        tree_score.canonicalize();
        return tree_score;
    }
}
