#!/usr/bin/env python
import GTsetMBparam, os, sys, argparse, subprocess

def main(argv):    
    # Set up variables for this program
    turnerdir = "Turner99"
    outputdir = "output/data"

    # Set up parameters
    parser = argparse.ArgumentParser(description="Calculate multibranch profile for an MFE structure")
    parser.add_argument("-s", "--structure", nargs=1, help="Structure to profile", required=True)

    args = vars(parser.parse_args())
    structfile = args["structure"][0]

    result = find_xyzw(turnerdir, outputdir, structfile)

def find_xyzw(turnerdir, outputdir, structfile):
    x = run_scorer(turnerdir, outputdir, structfile, [1, 0, 0, 0])
    y = run_scorer(turnerdir, outputdir, structfile, [0, 1, 0, 0])
    z = run_scorer(turnerdir, outputdir, structfile, [0, 0, 1, 0])
    w = run_scorer(turnerdir, outputdir, structfile, [0, 0, 0, 1])

    return [x, y, z, w]

def run_scorer(turnerdir, outputdir, seqfile, paramvec):
    # First, we set up an environment with the specified parameters
    GTsetMBparam.setup_gt_from_vec(turnerdir, outputdir, paramvec)

    # Then we run RNAScoring using those parameters
    result = subprocess.check_output(["RNAScoring", "--param-dir", os.path.split(outputdir)[0] + "/", seqfile])

    # The last line of the output contains the desired score
    lines = result.decode("utf-8").splitlines()
    lastline = lines[-1]
    score = round(float(lastline.split()[3]), 2)
    return score
    
# Voodoo to make Python run the program
if __name__ == "__main__":
    main(sys.argv[1:])

