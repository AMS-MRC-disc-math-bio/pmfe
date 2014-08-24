#!/usr/bin/env python
import os, sys, argparse, subprocess, shutil, logging, string
import RNAscorer
from gtmfe import gtmfe
from iB4e import iB4e
from parametrizer_types import parametrizer_types

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

    print run_iB4e(seqfile, sagefile, paramdir, structdir, verbose)

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
    def BBfunc(objectiveEV):
        param_vec = objectiveEV.as_param_vector()
        params_as_pairs = param_vec.get_pairs()
        pair_printer = lambda pair: str(pair[0]) + "/" + str(pair[1])
        
        result_file =  os.path.join(structdir, seqfilebase) + "." + "[" + " , ".join(pair_printer(pair) for pair in params_as_pairs) + "].ct"

        result_scores = gtmfe.mfe_main(seqfile, result_file, paramdir, param_vec)

        return iB4e.EuclideanVector(result_scores)
        
    return (classical_scores_python, BBfunc)

def run_iB4e(seqfile, sagefile, paramdir, structdir, verbose=False):
    logger = logging.getLogger()
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    if not (os.path.isfile(seqfile) and os.access(seqfile, os.R_OK)):
        raise IOError("Could not locate sequence file " + seqfile + ".")

    # Create the GTfold-running function
    (classical_scores_python, BBfunc) = setup_gtmfe_as_BlackBoxOptimize(seqfile, structdir, paramdir)

    # Set up a class to interface iB4e with GTfold
    class GTfoldPolytope(iB4e.BBPolytope):
        def BlackBoxOptimize(self, objective):
            return BBfunc(objective)

    thepoly = GTfoldPolytope(4)

    # Do magic
    thepoly.Build()

    # Retrieve the results
    vertices = thepoly.vertices

    # Now build a Sage file encoding the desired data
    build_sage_polytope_file(classical_scores_python, vertices, sagefile)

    # Finally, return the points for testing
    return vertices

def build_sage_polytope_file(classical_scores, vertices, sagefile):
    templatefile = open("output.template")
    template = string.Template(templatefile.read())

    point_pair_writer = lambda point_pair: "QQ(" + str(point_pair[0]) + "/" + str(point_pair[1]) + ")"
    pair_vector_writer = lambda pair_vector: "[" + " , ".join(point_pair_writer(pair) for pair in pair_vector) + "]"
    point_string = "[" + " , ".join(pair_vector_writer(vertex.as_param_vector().get_pairs()) for vertex in vertices) + "]"

    results = {"points": point_string, "classical_scores": classical_scores}

    sagecode = template.substitute(results)

    file = open(sagefile, mode='w+')
    file.write(sagecode)
    file.close()

    logging.info("Wrote Sage polytope file " + str(sagefile))

def score_parser(result):
    results_dict = result.get_python_fractions_dict()
    multiloops = results_dict["multiloops"]
    unpaired = results_dict["unpaired"]
    branches = results_dict["branches"]
    w = results_dict["w"]
    return [multiloops, unpaired, branches, w]

def find_w(scores, energy, params=[3.4, 0.0, 0.4, 1]):
    if params[3] == 0:
        return 0
    else:
        return (-scores[0]*params[0] + -scores[1]*params[1] + -scores[2]*params[2] + energy)/params[3]

# Voodoo to make Python run the program
if __name__ == "__main__":
    main(sys.argv[1:])
