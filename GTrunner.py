#!/usr/bin/env python
import GTsetMBparam, GTscorer, os, sys, argparse, subprocess

def main(argv):
    # Set up variables for this program
    turnerdir = "Turner99"
    outputdir = "output/data"
    scoredir = "output/scoring"
    
    # Set up parameters
    parser = argparse.ArgumentParser(description="Run GTfold with specified a, b, c, d parameters for iB4e")
    parser.add_argument("-s", "--sequence", nargs=1, help="Sequence to fold", required=True)
    parser.add_argument("-p", "--paramfile", nargs='?', help="Multibranch vector file (or give parameters on stdin", default=sys.stdin, type=argparse.FileType('r'))

    args = vars(parser.parse_args())
    paramfile = args["paramfile"]
    seqfile = args["sequence"][0]

    # First, we set up an environment for GTfold
    GTsetMBparam.setup_gt_from_file(turnerdir, outputdir, paramfile)

    structfile = run_gt(turnerdir, outputdir, seqfile)
    print "Structure stored in " + structfile

def run_gt(turnerdir, outputdir, seqfile):
    # Then we run GTfold on the specified sequence
    subprocess.check_output(["gtmfe", "-p", os.path.join(outputdir, turnerdir), seqfile])

    structfile = os.path.splitext(seqfile)[0] + ".ct"
    return structfile

# Voodoo to make Python run the program
if __name__ == "__main__":
    main(sys.argv[1:])

