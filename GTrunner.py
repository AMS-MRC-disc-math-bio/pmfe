#!/usr/bin/env python
import GTsetMBparam, os, sys, argparse, subprocess

def main(argv):
    # Set up variables for this program
    turnerdir = "Turner99"
    outputdir = "output"
    
    # Set up parameters
    parser = argparse.ArgumentParser(description="Run GTfold with specified a, b, c, d parameters for iB4e")
    parser.add_argument("-s", "--sequence", nargs=1, help="Sequence to fold", required=True)
    parser.add_argument("-p", "--paramfile", nargs='?', help="Multibranch vector file (or give parameters on stdin", default=sys.stdin, type=argparse.FileType('r'))

    args = vars(parser.parse_args())
    paramfile = args["paramfile"]
    seqfile = args["sequence"][0]

    run_gt(turnerdir, outputdir, paramfile, seqfile)

def run_gt(turnerdir, outputdir, paramfile, seqfile):
    # First, we set up an environment for GTfold
    GTsetMBparam.setup_gt(turnerdir, outputdir, paramfile)

    # Then we run GTfold on the specified sequence
    print outputdir
    subprocess.call(["gtmfe", "-p", os.path.join(outputdir, turnerdir), seqfile])

# Voodoo to make Python run the program
if __name__ == "__main__":
    main(sys.argv[1:])

