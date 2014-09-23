#!/usr/bin/env sage 
from sage.all import *
import os, sys, argparse, subprocess, shutil, logging, string
import RNAscorer
import BBpolytope
from gtmfe import gtmfe

def main(argv):
    # Set up parameters
    parser = argparse.ArgumentParser(description="Construct the branching polytope for a given RNA sequence")
    parser.add_argument("sequence", help="Sequence to fold")
    parser.add_argument("-v", "--verbose", help="Output debugging information", action="store_true")

    args = parser.parse_args()
    seqfile = args.sequence
    sagefile = os.path.splitext(seqfile)[0] + ".polytope.sage"
    verbose = args.verbose

    paramdir = "Turner99"

    structdir = os.path.splitext(seqfile)[0]
    try:
        os.makedirs(structdir)
    except OSError:
        if not os.path.isdir(structdir):
            raise

    run_iB4e(seqfile, sagefile, paramdir, structdir, verbose)

def run_iB4e(seqfile, sagefile, paramdir, structdir, verbose=False):
    logger = logging.getLogger()
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    if not (os.path.isfile(seqfile) and os.access(seqfile, os.R_OK)):
        raise IOError("Could not locate sequence file " + seqfile + ".")

    # Create the GTfold-running function
    (classical_scores, find_mfe_score) = setup_gtmfe_as_BlackBoxOptimize(seqfile, structdir, paramdir)

    # Do the magic
    thepoly = BBpolytope.build_polytope(find_mfe_score, 4)

    # Now build a Sage file encoding the desired data
    build_sage_polytope_file(classical_scores, thepoly, sagefile)

    # Finally, return the polytope for testing
    return thepoly

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
    def find_mfe_score(objective_vector):
        logging.debug("Optimizing for vector: " + str(objective_vector))
        
        parameter_vector = gtmfe.ParameterVector()
        pairify = lambda term: gtmfe.pairll(long(term.numerator()), long(term.denominator()))
        parameter_vector.set_from_pairs(pairify(objective_vector[0]), pairify(objective_vector[1]), pairify(objective_vector[2]), pairify(objective_vector[3]))

        result_file =  os.path.join(structdir, seqfilebase) + ".ct"
        result_scores = gtmfe.mfe_main(seqfile, result_file, paramdir, parameter_vector).get_python_fractions_dict()

        qqify = lambda frac: QQ(frac.numerator) / QQ(frac.denominator)
        result_vector = [qqify(result_scores["multiloops"]), qqify(result_scores["unpaired"]), qqify(result_scores["branches"]), qqify(result_scores["w"])]

        logging.debug("MFE structure scores: " + str(result_vector))

        result_vector_filename_string = "[ " + ", ".join([str(value.numerator()) for value in result_vector[:3]] + [str(result_vector[3].n(digits=3))]) + " ]"

        new_result_file = os.path.join(structdir, seqfilebase) + "." + result_vector_filename_string + ".ct"
        shutil.move(result_file, new_result_file)

        return result_vector

    return (classical_scores_gtmfe, find_mfe_score)

def build_sage_polytope_file(classical_scores, polytope, sagefile):
    templatefile = open("output.template")
    template = string.Template(templatefile.read())

    #point_pair_writer = lambda point_pair: "QQ(" + str(point_pair[0]) + "/" + str(point_pair[1]) + ")"
    #pair_vector_writer = lambda pair_vector: "[" + " , ".join(point_pair_writer(pair) for pair in pair_vector) + "]"
    #point_string = "[" + " , ".join(pair_vector_writer(vertex.as_param_vector().get_pairs()) for vertex in vertices) + "]"
    point_string = str(polytope.vertices_list())

    #classical_scores_pairs = [(classical_scores[0], 1), (classical_scores[1], 1), (classical_scores[2], 1), (classical_scores[3].numerator, classical_scores[3].denominator)]
    #classical_scores_string = pair_vector_writer(classical_scores_pairs)
    classical_scores_string = str(classical_scores)

    results = {"points": point_string, "classical_scores": classical_scores_string}

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
