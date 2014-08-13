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
    parser.add_argument("-a", help="Value of a (multibranch loop parameter)", type=float, default=3.4)
    parser.add_argument("-b", help="Value of b (unpaired nucleotide parameter)", type=float, default=0.0)
    parser.add_argument("-c", help="Value of c (branching helix parameter)", type=float, default=0.4)
    parser.add_argument("-d", help="Value of d (dummy scaling parameter)", type=float, default=1)

    args = parser.parse_args()
    seqfile = args.sequence
    verbose = args.verbose
    params = (args.a, args.b, args.c, args.d)

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

    print "a = {0}, b = {1}, c = {2}, d = {3}".format(params[0], params[1], params[2], params[3])
    print "x = {0}, y = {1}, z = {2}, w = {3}".format(result.multiloops, result.unpaired, result.branches, result.w)

def run_gtmfe(seqfile, structtarget, paramdir, params=[3.4, 0.0, 0.4, 1]):
    result = gtmfe.mfe_main(seqfile, structtarget, paramdir, params[0], params[1], params[2], params[3])
    return result

# Voodoo to make Python run the program
if __name__ == "__main__":
    main(sys.argv[1:])
