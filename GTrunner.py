#!/usr/bin/env python
import GTsetMBparam, GTscorer, os, sys, argparse, subprocess, logging

gtmfe_path = "gtmfe"

def main(argv):
    # Set up variables for this program
    inputdir = "input/gt/Turner99"
    outputdir = "output/gt/Turner99"
    scoredir = "output/scoring"
    
    # Set up parameters
    parser = argparse.ArgumentParser(description="Run GTfold with specified a, b, c, d parameters for iB4e")
    parser.add_argument("-s", "--sequence", nargs=1, help="Sequence to fold", required=True)
    parser.add_argument("-p", "--paramfile", nargs='?', help="Multibranch vector file (or give parameters on stdin", default=sys.stdin, type=argparse.FileType('r'))
    parser.add_argument("-v", "--verbose", help="Output debugging information", action="store_true")

    args = vars(parser.parse_args())
    paramfile = args["paramfile"]
    seqfile = args["sequence"][0]
    verbose = args["verbose"]

    logger = logging.getLogger()
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    # First, we set up an environment for GTfold
    
    GTsetMBparam.setup_gt_from_file(inputdir, outputdir, paramfile)

    structfile = run_gt(inputdir, outputdir, seqfile)
    print "Structure stored in " + structfile

def run_gt(inputdir, outputdir, seqfile=False):
    # Then we run GTfold on the specified sequence
    try:
        subprocess.check_output([gtmfe_path, "-p", outputdir, seqfile])
    except OSError:
        raise OSError("gtmfe executable not found!\nEdit the variable in " + __file__ + ".")

    logging.debug("Running gtmfe with path " + str(gtmfe_path) + " on sequence in " + str(seqfile))

    structfile = os.path.splitext(seqfile)[0] + ".ct"
    return structfile

# Voodoo to make Python run the program
if __name__ == "__main__":
    main(sys.argv[1:])

