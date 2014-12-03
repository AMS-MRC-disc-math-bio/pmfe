/*
 A class to score an RNA structure reading in from a file.
 */

#ifndef RNASCORING_H
#define RNASCORING_H

#include <string>
#include <gmpxx.h>
#include "StructureReader.h"
#include "RNAScoring.h"
#include "TreeScoring.h"

namespace rnascoring
{

    mpq_class get_classical_score(std::string structfile, std::string paramdir);

}

#endif
