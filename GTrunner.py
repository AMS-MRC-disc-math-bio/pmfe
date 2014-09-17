#!/usr/bin/env python
import os, sys, argparse, logging
from gtmfe import gtmfe

'''
Wrapper script to run gtmfe folding algorithm.
'''

def main(argv):
    # Set up parameters
    parser = argparse.ArgumentParser(description="Run GTfold with specified a, b, c, d parameters")
    parser.add_argument("sequence", help="Sequence to fold")
    parser.add_argument("-v", "--verbose", help="Output debugging information", action="store_true")
    parser.add_argument("-o", "--structure", nargs=1, help="File to store structure result")
    parser.add_argument("-p", "--parameters", nargs=4, help="List of parameters", metavar=("A", "B", "C", "D"), default=("3.4", "0.0", "0.4", "1"))

    args = parser.parse_args()
    seqfile = args.sequence
    verbose = args.verbose
    params = gtmfe.ParameterVector()
    params.set_from_words(args.parameters[0], args.parameters[1], args.parameters[2], args.parameters[3])

    logger = logging.getLogger()
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
        
    paramdir = "Turner99"

    # Use the supplied output file if applicable
    try:
        structtarget = args.structure[0]
        # Otherwise, generate one from the input filename
    except TypeError:
        structtarget = os.path.splitext(os.path.basename(seqfile))[0] + ".ct"

    result = run_gtmfe(seqfile, structtarget, paramdir, params)

    print "Results: " + str(score_parser(result))

def run_gtmfe(seqfile, structtarget, paramdir, params=gtmfe.ParameterVector()):
    result = gtmfe.mfe_main(seqfile, structtarget, paramdir, params)
    return result

def score_parser(result):
    results_dict = result.get_python_fractions_dict()
    multiloops = results_dict["multiloops"]
    unpaired = results_dict["unpaired"]
    branches = results_dict["branches"]
    w = results_dict["w"]
    return [multiloops, unpaired, branches, w]

# Voodoo to make Python run the program
if __name__ == "__main__":
    main(sys.argv[1:])
