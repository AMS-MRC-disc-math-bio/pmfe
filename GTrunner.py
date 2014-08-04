#!/usr/bin/env python
import os, sys, argparse, logging
from gtmfe import gtmfe

def main(argv):
    # Set up parameters
    parser = argparse.ArgumentParser(description="Run GTfold with specified a, b, c, d parameters")
    parser.add_argument("-s", "--sequence", nargs=1, help="Sequence to fold", required=True)
    parser.add_argument("-v", "--verbose", help="Output debugging information", action="store_true")
    parser.add_argument("-a", help="Value of a", type=float, default=10.1)
    parser.add_argument("-b", help="Value of b", type=float, default=0.3)
    parser.add_argument("-c", help="Value of c", type=float, default=0.3)
    parser.add_argument("-d", help="Value of d", type=float, default=1)

    args = vars(parser.parse_args())
    seqfile = args["sequence"][0]
    verbose = args["verbose"]
    params = (args["a"], args["b"], args["c"], args["d"])

    logger = logging.getLogger()
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
        
    paramdir = "Turner99"

    structdir = os.path.splitext(seqfile)[0]
    structname = os.path.splitext(os.path.basename(seqfile))[0] + "." + str(params) + ".ct"
    structtarget = os.path.join(structdir, structname)
    result = gtmfe.mfe_main(seqfile, structtarget, paramdir, params[0], params[1], params[2], params[3])

    print "a = {0}, b = {1}, c = {2}, d = {3}".format(params[0], params[1], params[2], params[3])
    print "x = {0}, y = {1}, z = {2}, w = {3}".format(result.x, result.y, result.z, result.w)
    

# Voodoo to make Python run the program
if __name__ == "__main__":
    main(sys.argv[1:])
