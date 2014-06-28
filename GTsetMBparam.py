#!/usr/bin/env python
import sys, argparse, os, shutil, fileinput, re, numpy, logging

def main(argv):
    # Set up and process arguments   
    parser = argparse.ArgumentParser(description="Set up GTfold with specified multibranch energy parameters")
    parser.add_argument("-t", "--turnerdir", nargs=1, help="Location of original Turner99 DAT files", required=True)
    parser.add_argument("-o", "--outputdir", nargs=1, help="Output base directory to receive new Turner99 directory", required=True)
    parser.add_argument("-p", "--paramfile", nargs='?', help="Multibranch vector file (or give parameters on stdin", default=sys.stdin, type=argparse.FileType('r'))
    parser.add_argument("-v", "--verbose", help="Output debugging information", action="store_true")
    
    args = vars(parser.parse_args())
    paramfile = args["paramfile"]
   
    turnerdir = os.path.abspath(args["turnerdir"][0])
    outputdir = os.path.abspath(args["outputdir"][0])
    verbose = args["verbose"]

    logger = logging.getLogger()
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    setup_gt_from_file(turnerdir, outputdir, paramfile)
    
def setup_gt_from_file(turnerdir, outputdir, paramfile):
    # Read the numerical values of the parameters
    new_params = get_params(paramfile)

    setup_gt_from_vec(turnerdir, outputdir, new_params)

def copy_turner(turnerdir, targetdir):
    logging.debug("Copied Turner99 parameters to " + str(targetdir))
    shutil.rmtree(targetdir, ignore_errors=True)
    shutil.copytree(turnerdir, targetdir)

def setup_gt_from_vec(turnerdir, outputdir, new_params):
    # We separate out the program logic so this can be run from other Python scripts
    # Define the target directory name
    targetdir = os.path.join(outputdir, "Turner99")

    # Copy the Turner99 data files into the new directory
    copy_turner(turnerdir, targetdir)
     
    # Modify the copied Turner99 files with the parameters
    write_new_params(targetdir, new_params)
    logging.debug("Wrote parameters " + str(new_params) + " to " + str(targetdir))
    
def get_params(paramfile):
    # Note that we negate the parameters since GTfold is a minimizer
    params = [-round(float(p), 2) for p in paramfile.readline().split()]
    return params

def write_new_params(targetdir, new_params):
    # Locate the Turner99 files and build a generator of their absolute paths
    files = (os.path.join(targetdir, file) for file in  os.listdir(targetdir))
    miscfile = os.path.join(targetdir, "Turner99")

    normalized_params = normalize_params(new_params)

    # Open all the files and munge them using the new parameters
    input = fileinput.input(files, inplace=True)
    for line in input:
        rewrite_line(line, normalized_params)
    input.close()

def normalize_params(parameters):
    p = numpy.float_(parameters)
    pnorm = numpy.linalg.norm(p)
    
    if pnorm > 100:
        result = numpy.rint(100*p/pnorm)
    else:
        result = p
    return result

def rewrite_line(line, new_params):
    # Each line of the file is made of up "words", which may be numbers or other string data.
    # We first split them into a list.
    words = re.split(r'(\s+)', line.rstrip('\n'))

    # If this is the multibranch parameters, replace them with a, b, c
    if "3.40" in words and ".00" in words and ".40" in words:
        newwords = [process_word_abc(word, new_params[:3]) for word in words]
    # Otherwise, multiply all the numerical words by d
    else:
        newwords = [process_word_generic(word, new_params[3]) for word in words]

    # Finally, we build a new line separated by tab characters.
    # The semantics of fileinput require that we print the desired line to stdout.
    print ''.join(newwords)

def process_word_generic(word, d):
    try:
        # Leave the word alone if it's an integer
        newword = str(int(word))

    except ValueError:
        try:
            if word == "inf":
                # Leave inf as it is, no matter what
                newword = word
            else:
                # Multiply the word by d if it's a float
                newword = "{:.3f}".format(d*float(word))[:len(word)]

        except ValueError:
            # Some of the words are strings of nucleotides, which can't be cast to float
            newword = word
            pass
            
    return newword

def process_word_abc(word, params):
    a, b, c = params
    if word == "3.40":
        newword = '%.2f' % float(a) 
    elif word == ".00":
        newword = '%.2f' % float(b) 
    elif word == ".40":
        newword = '%.2f' % float(c) 
    else:
        newword = word

    return str(newword)

# Voodoo to make Python run the program
if __name__ == "__main__":
    main(sys.argv[1:])
