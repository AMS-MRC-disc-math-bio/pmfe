#!/usr/bin/env python
import GTsetMBparam, GTscorer, os, sys, argparse, subprocess

def main(argv):
    # Set up variables for this program
    turnerdir = "Turner99"
    outputdir = "output/data"
    
    # Set up parameters
    parser = argparse.ArgumentParser(description="Run GTfold with specified a, b, c, d parameters for iB4e")
    parser.add_argument("-s", "--sequence", nargs=1, help="Sequence to fold", required=True)
    parser.add_argument("-o", "--polyfile", nargs=1, help="Name for polymake output", required=True)
    parser.add_argument("-p", "--paramfile", nargs='?', help="Multibranch vector file (or give parameters on stdin", default=sys.stdin, type=argparse.FileType('r'))

    args = vars(parser.parse_args())
    paramfile = args["paramfile"]
    polyfile = args["polyfile"][0]
    seqfile = args["sequence"][0]

    # First, we set up an environment for GTfold
    GTsetMBparam.setup_gt_from_file(turnerdir, outputdir, paramfile)

    run_gt(turnerdir, outputdir, seqfile)

def run_gt(turnerdir, outputdir, seqfile):
    # Then we run GTfold on the specified sequence
    subprocess.check_output(["gtmfe", "-p", os.path.join(outputdir, turnerdir), seqfile])

    structfile = os.path.splitext(seqfile)[0] + ".ct"
    return structfile
    
def run_scorer(turnerdir, outputdir, structfile):
    x, y, z, w = GTscorer.find_xyzw(turnerdir, outputdir, structfile)
    return [x, y, z, w]

def write_polyfile_vector(vector, polyfile):
    # Write out the vector in Polymake format with homogenous coordinates
    file = open(polyfile, mode='a')
    vector_string = "[1," + ",".join(map(str, vector)) + "],"
    file.write(vector_string)
    file.close()

# Voodoo to make Python run the program
if __name__ == "__main__":
    main(sys.argv[1:])

