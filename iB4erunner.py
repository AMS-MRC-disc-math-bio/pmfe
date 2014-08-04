#!/usr/bin/env python
import os, sys, argparse, subprocess, shutil, logging, string
import RNAscorer
from gtmfe import gtmfe

iB4e_path = "iB4e/iB4e-rna"

def main(argv):    
    # Set up parameters
    parser = argparse.ArgumentParser(description="Run iB4e and GTFold to construct the structure polytope")
    parser.add_argument("-s", "--sequence", nargs=1, help="Sequence to fold", required=True)
    parser.add_argument("-v", "--verbose", help="Output debugging information", action="store_true")

    args = vars(parser.parse_args())
    seqfile = args["sequence"][0]
    sagefile = os.path.splitext(seqfile)[0] + ".sage"
    verbose = args["verbose"]

    logger = logging.getLogger()
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
        

    paramdir = "Turner99"
    
    structdir = os.path.splitext(seqfile)[0]
    try:
        os.makedirs(structdir)
    except OSError:
        if not os.path.isdir(structdir):
            raise
            
    points = []

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
    classical_file = os.path.join(structdir, os.path.splitext(seqfile)[0] + ".classical.ct")
    classical_result = score_parser(gtmfe.mfe_main(seqfile, classical_file, paramdir))
    classical_scores = RNAscorer.score_file(classical_file)

    logging.debug("Classical scores from Python: " + str(classical_scores))
    logging.debug("Classical scores from gtmfe: " + str(classical_result))

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
        structname = os.path.splitext(os.path.basename(seqfile))[0] + "." + str(params) + ".ct"
        structtarget = os.path.join(structdir, structname)
        result = score_parser(gtmfe.mfe_main(seqfile, structtarget, paramdir, params[0], params[1], params[2], params[3]))
        logging.debug("Stored structure as " + str(structtarget))
        
        # Store the scores
        scores = RNAscorer.score_file(structtarget)
        if scores != (result[0], result[1], result[2]):
            logging.warn("Score mismatch for parameters {0}!".format(params))
            logging.warn("GTfold gives x={0}, y={1}, z={2}.".format(result[0], result[1], result[2]))
            logging.warn("RNAscorer gives x={0}, y={1}, z={2}".format(scores[0], scores[1], scores[2]))
            
        points.append(scores)
        result = " ".join(map(str, scores)) + "\n"
        logging.debug("Structure scores " + str(scores))

        # Send the score vector to iB4e
        iB4e.stdin.write(result)

    # Now build a Sage file encoding the desired data
    build_sage_polytope_file(classical_scores, points, sagefile)
    
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

# Voodoo to make Python run the program
if __name__ == "__main__":
    main(sys.argv[1:])
