#!/usr/bin/env python
import os, sys, argparse, subprocess, shutil, logging, string
import RNAscorer
from gtmfe import gtmfe

iB4e_path = "iB4e/iB4e-rna"

def main(argv):
    # Set up parameters
    parser = argparse.ArgumentParser(description="Run iB4e and GTFold to construct the structure polytope")
    parser.add_argument("sequence", help="Sequence to fold")
    parser.add_argument("-v", "--verbose", help="Output debugging information", action="store_true")

    args = parser.parse_args()
    seqfile = args.sequence
    sagefile = os.path.splitext(seqfile)[0] + ".sage"
    verbose = args.verbose

    paramdir = "Turner99"

    structdir = os.path.splitext(seqfile)[0]
    try:
        os.makedirs(structdir)
    except OSError:
        if not os.path.isdir(structdir):
            raise

    run_iB4e(seqfile, sagefile, paramdir, structdir, verbose=False)

def run_iB4e(seqfile, sagefile, paramdir, structdir, verbose=False):
    logger = logging.getLogger()
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    if not (os.path.isfile(seqfile) and os.access(seqfile, os.R_OK)):
        raise IOError("Could not locate sequence file " + seqfile + ".")

    points = []

    seqfilebase = os.path.splitext(os.path.basename(seqfile))[0]

    logging.debug("Starting iB4e")

    try:
        iB4e = subprocess.Popen([iB4e_path, '4'], # Our problem operates in 4 dimensions
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                shell=False)
    except OSError:
        raise OSError("iB4e-rna executable not found!\nBe sure it is built, then edit the variable in " + __file__ + ".")

    logging.debug("Computing classical scores")

    # Compute classical scores for later reference
    classical_file = os.path.join(structdir, seqfilebase + ".classical.ct")
    classical_result = gtmfe.mfe_main(seqfile, classical_file, paramdir)
    classical_scores_gtmfe = score_parser(classical_result)
    classical_scores_python = list(RNAscorer.score_file(classical_file))
    
    classical_scores_python.append(find_w(classical_scores_python, classical_result.energy))

    logging.debug("Classical scores from Python: " + str(classical_scores_python))
    logging.debug("Classical scores from gtmfe: " + str(classical_scores_gtmfe))

    # Create storage directory for structures
    shutil.rmtree(structdir, ignore_errors=True)
    os.mkdir(structdir)

    # Now turn iB4e loose to compute the geometry
    while True:
        # At the end of the file, we'll receive a blank line, which makes eval() choke
        try:
            params = eval(iB4e.stdout.readline())
        except (SyntaxError, TypeError):
            break

        logging.debug("iB4e requests vector " + str(params))

        # iB4e is a maximizer, while gtmfe is a minimizer
        params = tuple([-p for p in params])

        # Find the MFE structure
        result_file =  os.path.join(structdir, seqfilebase) + "." + str(params) + ".ct"
        result_gtmfe = gtmfe.mfe_main(seqfile, result_file, paramdir, params[0], params[1], params[2], params[3])
        scores_gtmfe = score_parser(result_gtmfe)
        logging.debug("Stored structure as " + str(result_file))

        # Store the scores
        scores_python = list(RNAscorer.score_file(result_file))
        if scores_python != scores_gtmfe[:3]:
            logging.warn("Score mismatch for parameters {0}!".format(params))
            logging.warn("GTfold gives x={0}, y={1}, z={2}.".format(scores_gtmfe[0], scores_gtmfe[1], scores_gtmfe[2]))
            logging.warn("RNAscorer gives x={0}, y={1}, z={2}".format(scores_python[0], scores_python[1], scores_python[2]))
        scores_python.append(find_w(scores_python, result_gtmfe.energy, params))

        points.append(scores_python)
        result = " ".join(map(str, scores_python)) + "\n"
        logging.debug("Structure scores from Python: " + str(scores_python))
        logging.debug("Structure scores from gtmfe: " + str(scores_gtmfe))

        # Send the score vector to iB4e
        iB4e.stdin.write(result)

    # Now build a Sage file encoding the desired data
    build_sage_polytope_file(classical_scores_python, points, sagefile)

    # Finally, return the points for testing
    return points

def build_sage_polytope_file(classical_scores, points, sagefile):
    templatefile = open("output.template")
    template = string.Template(templatefile.read())

    results = {"points": points, "classical_scores": classical_scores}

    sagecode = template.substitute(results)

    file = open(sagefile, mode='w+')
    file.write(sagecode)
    file.close()

    logging.info("Wrote Sage polytope file " + str(sagefile))

def score_parser(result):
    return [int(result.multiloops), int(result.unpaired), int(result.branches), result.w]

def find_w(scores, energy, params=[3.4, 0.4, 0.0, 1]):
    if params[3] == 0:
        return 0
    else:
        return (-scores[0]*params[0] + -scores[1]*params[1] + -scores[2]*params[2] + energy)/params[3]

# Voodoo to make Python run the program
if __name__ == "__main__":
    main(sys.argv[1:])
