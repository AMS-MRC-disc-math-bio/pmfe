#!/usr/bin/env python
import GTsetMBparam, os, sys, argparse, subprocess

def main(argv):    
    # Set up variables for this program
    turnerdir = "Turner99"
    scoredir = "output/scoring"

    # Set up parameters
    parser = argparse.ArgumentParser(description="Calculate multibranch profile for an MFE structure")
    parser.add_argument("-s", "--structure", nargs=1, help="Structure to profile", required=True)

    args = vars(parser.parse_args())
    structfile = args["structure"][0]

    result = find_xyzw(turnerdir, scoredir, structfile)
    print result

def find_xyzw(turnerdir, scoredir, structfile):
    x = run_vscorer(turnerdir, scoredir, "x", [1,0,0,0], structfile, as_float=False)
    y = run_vscorer(turnerdir, scoredir, "y", [0,1,0,0], structfile, as_float=False)
    z = run_vscorer(turnerdir, scoredir, "z", [0,0,1,0], structfile, as_float=False)
    w = run_vscorer(turnerdir, scoredir, "w", [0,0,0,1], structfile, as_float=True)

    return [x, y, z, w]

def run_vscorer(turnerdir, scoredir, vname, paramvec, structfile, as_float=True):
    vdir = os.path.join(scoredir, vname, "data")
    if not os.path.isfile(os.path.join(vdir, "Turner99", "miscloop.dat")):
        setup_scorer(turnerdir, vdir, paramvec)

    return run_scorer(vdir, structfile, as_float)

def setup_scorer(turnerdir, outputdir, paramvec):
    # First, we set up an environment with the specified parameters
    GTsetMBparam.setup_gt_from_vec(turnerdir, outputdir, paramvec)

def run_scorer(outputdir, structfile, as_float=True):
    result = subprocess.check_output(["RNAScoring", "--param-dir", os.path.split(outputdir)[0] + "/", structfile])

    # The last line of the output contains the desired score
    lines = result.decode("utf-8").splitlines()
    lastline = lines[-1]
    if as_float:
        score = round(float(lastline.split()[3]), 2)
    else:
        # Due to a glitch in RNAScoring, we sometimes get spurious noise below the decimal point.
        # Under the assumption that there are never more than 99 poly-C hairpin loops,
        # we just drop this information.
        score = int(float(lastline.split()[3]))
    return score
    
# Voodoo to make Python run the program
if __name__ == "__main__":
    main(sys.argv[1:])

