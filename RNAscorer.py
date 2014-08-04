#!/usr/bin/env python

import os, sys, argparse, logging
from tree import PartitionIntervalNode

def main(argv):
    # Set up parameters
    parser = argparse.ArgumentParser(description="Score a specified RNA secondary structure")
    parser.add_argument("-c", "--structure", nargs=1, help="Structure to score (in CT format)", required=True)
    parser.add_argument("-v", "--verbose", help="Output debugging information", action="store_true")

    args = vars(parser.parse_args())
    structfile = args["structure"][0]
    verbose = args["verbose"]

    logger = logging.getLogger()
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    score_file(structfile)

def score_file(structfile):
    pairs, seqlength = parse_ct_file(structfile)
    tree = build_tree(pairs, seqlength)
    scores = score_tree(tree, is_root=True)
    logging.debug("Tree scores: x = {0}, y = {1}, z = {2}".format(scores[0], scores[1], scores[2]))
    return scores

# Read in a structure in CT format and produce a list of ordered pairs
def parse_ct_file(structfile):
    logging.debug("Opening structure file " + str(structfile))
    f = open(structfile, 'r')
    
    def lineparser(line):
        # Lines not representing bases are short, so we can catch them using IndexError
        try:
            terms = line.split()
            i, j = int(terms[0]), int(terms[4])
        except IndexError:
            return None

        # Drop pairs involving 0 or which are out of order
        if i < j and j != 0:
            return (i, j)
        else:
            return None

    # The first line contains the number of pairs
    seqlength = int(f.readline().split()[0])
    
    # The rest of the lines encode pairings
    pairs = [lineparser(line) for line in f if lineparser(line) is not None]
    f.close()

    logging.debug("Parsed structure file successfully.")
    
    return pairs, seqlength

def build_tree(pairs, seqlength):
    tree = PartitionIntervalNode(0, seqlength+1)
    for pair in pairs:
        tree.insert(pair[0], pair[1])

    logging.debug("Sorting tree")
    tree.sort()
    return tree

def score_tree(tree, is_root=False):
    x = 0 # Multiloop counter
    y = 0 # Unpaired base counter
    z = 0 # Branched helix counter
    if tree.valency() >= 2 and not is_root:
        x = 1
        y = tree.end - tree.start - 1 - sum(child.end - child.start + 1 for child in tree.children)
        z = tree.valency() + 1

    child_scores = [score_tree(child) for child in tree.children]

    x += sum(scores[0] for scores in child_scores)
    y += sum(scores[1] for scores in child_scores)
    z += sum(scores[2] for scores in child_scores)

    return x, y, z

# Voodoo to make Python run the program
if __name__ == "__main__":
    main(sys.argv[1:])
