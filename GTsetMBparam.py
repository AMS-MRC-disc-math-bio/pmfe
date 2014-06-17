#!/usr/bin/env python

import sys, argparse, os, shutil

def main(argv):
    paramfile = ''
    turnerdir = ''
    outputdir = ''
    
    parser = argparse.ArgumentParser(description="Set up GTfold with specified multibranch energy parameters")

    parser.add_argument("-p", "--paramfile", nargs=1, help="Multibranch vector file location", required=True)
    parser.add_argument("-t", "--turnerdir", nargs=1, help="Location of original Turner99 DAT files", required=True)
    parser.add_argument("-o", "--outputdir", nargs=1, help="Output base directory to receive new Turner99 directory", required=True)
    
    args = vars(parser.parse_args())
    paramfile = os.path.abspath(args["paramfile"][0])
    turnerdir = os.path.abspath(args["turnerdir"][0])
    outputdir = os.path.abspath(args["outputdir"][0])
    targetdir = os.path.join(outputdir, "Turner99")

    copy_turner(turnerdir, targetdir)

    new_params = get_params(paramfile)

    write_params(targetdir, new_params)

def copy_turner(turnerdir, targetdir):
    shutil.rmtree(targetdir, ignore_errors=True)
    shutil.copytree(turnerdir, targetdir)
    
def get_params(paramfile):
    f = open(paramfile, 'r')
    params = [str(round(float(p), 2)) for p in f.readline().split()]
    return params

def write_params(targetdir, new_params):
    miscfile = os.path.join(targetdir, "miscloop.DAT")
    miscbak = os.path.join(targetdir, "miscloop.DAT.BAK")

    shutil.move(miscfile, miscbak)

    newline = "\t".join(new_params) + "\n"
    
    source = open(miscbak, 'r')
    target = open(miscfile, 'w')
    
    lines = source.readlines()
    lines[3] = newline # We believe the third line of this file controls the multibranch loops

    target.writelines(lines)    

if __name__ == "__main__":
    main(sys.argv[1:])
