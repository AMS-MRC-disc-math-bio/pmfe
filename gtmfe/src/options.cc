#include "loader.h"
#include "options.h"

#include <sstream>

using namespace std;

bool ILSA;
bool NOISOLATE;
bool PARAM_DIR = false;
bool LIMIT_DISTANCE;
bool BPP_ENABLED;
bool SUBOPT_ENABLED;
bool CONS_ENABLED = false;
bool VERBOSE = false;
//bool SHAPE_ENABLED = false;
bool T_MISMATCH = false;
bool UNAMODE = false;
bool RNAMODE = false;
bool b_prefilter = false;
bool CALC_PART_FUNC = false;
bool RND_SAMPLE = false;
bool PF_COUNT_MODE = false;

string seqfile = "";
string constraintsFile = "";
string outputPrefix = "";
string outputFile = "";
string suboptFile = "";
string bppOutFile = "";
string outputDir = "";
string shapeFile = "";
string paramDir; // default value

int num_rnd = 0;
int dangles=-1;
int prefilter1=2;
int prefilter2=2;

float suboptDelta = 0.0;
int nThreads = -1;
int contactDistance = -1;
