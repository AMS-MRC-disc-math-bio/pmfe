import pyparam
from fractions import Fraction
import itertools

# This example shows how we can run gtmfe iteratively on a large collection of scores

# For this example, we'll use the "test_tRNA.fasta" sequence in the "test_data" folder
seqfile = "test_data/test_tRNA.fasta"

# and write the result to "test.ct" in the base directory
outfile = "test.ct"

mfe_scores = []

# We use itertools.product to generate all 4-vectors with coordinates in [-1, 0, 1]
for vector in itertools.product([-1, 0, 1], repeat=4):
    # For each of these vectors, we generate the corresponding ParameterVector
    pv = pyparam.ParameterVector.from_py_params(vector)
    # and then call mfe_scores with it
    mfe_scores.append(pyparam.get_mfe_scores(seqfile, outfile, pv).as_fractions())
    ### The .as_fractions() method gets the result as a list of instances of Python's Fraction, a numeric type. A .as_rational() method is also available which gives instances of Sage's Rational class instead.

# In this case, printing all 81 results isn't super informative, but what the heck?
# In general, you'll probably have a more useful idea of what to do with these results (such as checking whether they lie inside a pre-generated polytope)
print "All done!"
print mfe_scores
