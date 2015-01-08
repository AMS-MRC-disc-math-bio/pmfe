import pyparam
from fractions import Fraction

# This example shows how to run gtmfe with a single specified set of parameters

# First, we construct a ParameterVector representing our parameters.
# This may include any combination of Python ints, Python Fractions, and Sage Rationals (as long as we're in Sage, of course!)
pv = pyparam.ParameterVector.from_py_params([1, -1, Fraction(3, 2), -1])

# For this example, we'll use the "test_tRNA.fasta" sequence in the "test_data" folder
seqfile = "test_data/test_tRNA.fasta"

# and write the result to "test.ct" in the base directory
outfile = "test.ct"

# Next, we get the MFE scores for this sequence under these parameters
mfe_scores = pyparam.get_mfe_scores(seqfile, outfile, pv)

# and print the results
print mfe_scores
