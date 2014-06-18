#!/usr/bin/env python
import sys, argparse, os, shutil, fileinput

def main(argv):
    # Set up and process arguments   
    parser = argparse.ArgumentParser(description="Set up GTfold with specified multibranch energy parameters")
    parser.add_argument("-t", "--turnerdir", nargs=1, help="Location of original Turner99 DAT files", required=True)
    parser.add_argument("-o", "--outputdir", nargs=1, help="Output base directory to receive new Turner99 directory", required=True)
    parser.add_argument("-p", "--paramfile", nargs='?', help="Multibranch vector file (or give parameters on stdin", default=sys.stdin, type=argparse.FileType('r'))
    
    args = vars(parser.parse_args())
    paramfile = args["paramfile"]
    turnerdir = os.path.abspath(args["turnerdir"][0])
    outputdir = os.path.abspath(args["outputdir"][0])
    targetdir = os.path.join(outputdir, "Turner99")

    # Copy the Turner99 data files into the new directory
    copy_turner(turnerdir, targetdir)

    # Make a vector of floats of the new parameters
    new_params = get_params(paramfile)

    # Modify the copied Turner99 files with the parameters
    write_new_params(targetdir, new_params)
        
def copy_turner(turnerdir, targetdir):
    shutil.rmtree(targetdir, ignore_errors=True)
    shutil.copytree(turnerdir, targetdir)
    
def get_params(paramfile):
    # Note that we negate the parameters since GTfold is a minimizer
    params = [-round(float(p), 2) for p in paramfile.readline().split()]
    return params

def write_new_params(targetdir, new_params):
    # Locate the Turner99 files and build a generator of their absolute paths
    files = (os.path.join(targetdir, file) for file in  os.listdir(targetdir))
    miscfile = os.path.join(targetdir, "Turner99")

    # Open all the files and munge them using the new parameters
    input = fileinput.input(files, inplace=True)
    for line in input:
        rewrite_line(line, new_params)
    input.close()

def rewrite_line(line, new_params):
    # Each line of the file is made of up "words", which may be numbers or other string data.
    # We first split them into a list.
    words = line.split()

    # If this is the multibranch parameters, replace them with a, b, c
    if words == ["3.40", ".00", ".40"]:
        newwords = map(str, new_params[:3])
    # Otherwise, multiply all the numerical words by d
    else:
        newwords = [process_word(word, new_params[3]) for word in words]

    # Finally, we build a new line separated by tab characters.
    # The semantics of fileinput require that we print the desired line to stdout.
    print "\t".join(newwords)

def process_word(word, d):
    try:
        # Multiply the word by d if it's a number
        word = d*float(word) 
    except ValueError:
        # Some of the words are strings of nucleotides, which can't be cast to float
        pass
    return str(word)

# Voodoo to make Python run the program
if __name__ == "__main__":
    main(sys.argv[1:])
