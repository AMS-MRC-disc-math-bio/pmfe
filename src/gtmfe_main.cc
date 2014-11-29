// Copyright (c) 2014 Andrew Gainer-Dewar.

#include "mfe.h"
#include "parametrizer_types.h"
#include <iostream>

int main() {
    std::cout << mfe("combinatorial_seq.fasta", "combinatorial_seq.ct", "Turner99") << std::endl;
    return(0);
}
