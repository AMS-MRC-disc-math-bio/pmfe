#!/usr/bin/env python
import os, sys, argparse, subprocess, shutil, logging, string
import RNAscorer
from gtmfe import gtmfe
from iB4e import iB4e

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

    print run_iB4e(seqfile, sagefile, paramdir, structdir, verbose=False)

def setup_gtmfe_as_BlackBoxOptimize(seqfile, structdir, paramdir):
    # Create storage directory for structures
    shutil.rmtree(structdir, ignore_errors=True)
    os.mkdir(structdir)
    seqfilebase = os.path.splitext(os.path.basename(seqfile))[0]

    logging.debug("Computing classical scores")

    # Compute classical scores for later reference
    classical_file = os.path.join(structdir, seqfilebase + ".classical.ct")
    classical_result = gtmfe.mfe_main(seqfile, classical_file, paramdir)
    classical_scores_gtmfe = score_parser(classical_result)
    classical_scores_python = list(RNAscorer.score_file(classical_file))
    
    #classical_scores_python.append(find_w(classical_scores_python, classical_result.energy))

    logging.debug("Classical scores from Python: " + str(classical_scores_python))
    logging.debug("Classical scores from gtmfe: " + str(classical_scores_gtmfe))

    # Build the function which wraps gtmfe input and output in iB4e-compatible structures
    def BBfunc(params_as_EV):
        params_as_pairs = params_as_EV.get_split_values()

        pair_printer = lambda pair: str(pair[0]) + "/" + str(pair[1])
        
        result_file =  os.path.join(structdir, seqfilebase) + "." + "[" + " , ".join(pair_printer(pair) for pair in params_as_pairs) + "].ct"

        pairs_as_pairll = [gtmfe.pairll(pair[0], pair[1]) for pair in params_as_pairs]
        params_as_PV = gtmfe.ParameterVector()
        params_as_PV.set_from_pairs(pairs_as_pairll[0], pairs_as_pairll[1], pairs_as_pairll[2], pairs_as_pairll[3])
        
        result_gtmfe = gtmfe.mfe_main(seqfile, result_file, paramdir, params_as_PV)
        result_pairs = result_gtmfe.get_pairs()

        result = iB4e.EuclideanVector(4)
        result.set_split_values(result_pairs)

        return result

    return BBfunc        

def run_iB4e(seqfile, sagefile, paramdir, structdir, verbose=False):
    logger = logging.getLogger()
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    if not (os.path.isfile(seqfile) and os.access(seqfile, os.R_OK)):
        raise IOError("Could not locate sequence file " + seqfile + ".")

    # Create the GTfold-running function
    BBfunc = setup_gtmfe_as_BlackBoxOptimize(seqfile, structdir, paramdir)

    # Set up a class to interface iB4e with GTfold
    class GTfoldPolytope(iB4e.BBPolytope):
        def BlackBoxOptimize(self, objective):
            return BBfunc(objective)

    thepoly = GTfoldPolytope(4)

    # Do magic
    thepoly.Build()

    # Retrieve the results
    point_pairs = thepoly.get_split_vertices()

    # Now build a Sage file encoding the desired data
    #build_sage_polytope_file(classical_scores_python, point_pairs, sagefile)

    # Finally, return the points for testing
    return point_pairs

def build_sage_polytope_file(classical_scores, point_pairs, sagefile):
    templatefile = open("output.template")
    template = string.Template(templatefile.read())

    point_pair_writer = lambda point_pair: "QQ(" + str(pair[0]) + "/" + str(pair[1]) + ")"
    pair_vector_writer = lambda pair_vector: "[" + " , ".join(point_pair_writer(pair) for pair in pair_vector) + "]"
    point_string = "[" + " , ".join(pair_vector_writer(pair_vector) for pair_vector in point_pairs) + "]"

    results = {"points": point_string, "classical_scores": classical_scores}

    sagecode = template.substitute(results)

    file = open(sagefile, mode='w+')
    file.write(sagecode)
    file.close()

    logging.info("Wrote Sage polytope file " + str(sagefile))

def score_parser(result):
    return [int(result.multiloops), int(result.unpaired), int(result.branches), result.w]

def find_w(scores, energy, params=[3.4, 0.0, 0.4, 1]):
    if params[3] == 0:
        return 0
    else:
        return (-scores[0]*params[0] + -scores[1]*params[1] + -scores[2]*params[2] + energy)/params[3]

# Voodoo to make Python run the program
if __name__ == "__main__":
    main(sys.argv[1:])
