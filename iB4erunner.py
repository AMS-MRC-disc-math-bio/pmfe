#!/usr/bin/env python
import os, sys, argparse, subprocess
import GTrunner, GTsetMBparam, GTscorer

def main(argv):    
    # Set up parameters
    parser = argparse.ArgumentParser(description="Run GTfold with specified a, b, c, d parameters for iB4e")
    parser.add_argument("-s", "--sequence", nargs=1, help="Sequence to fold", required=True)

    args = vars(parser.parse_args())
    seqfile = args["sequence"][0]
    sagefile = os.path.splitext(seqfile)[0] + ".sage"

    turnerdir = "Turner99"
    outputdir = "output/data"

    points = []

    iB4e = subprocess.Popen(['iB4e/iB4e', '4'],
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            shell=False)

    while True:
        # At the end of the file, we'll receive an EOF character, which makes eval() choke
        try:
            params = eval(iB4e.stdout.readline())
        except SyntaxError:
            break

        # Set up the parameters
        GTsetMBparam.setup_gt_from_vec(turnerdir, outputdir, params)

        # Find the MFE structure
        structfile = GTrunner.run_gt(turnerdir, outputdir, seqfile)
        
        # Score the structure using GTfold
        scores = GTscorer.find_xyzw(turnerdir, outputdir, structfile)
        points.append(scores)
        result = " ".join(map(str, scores)) + "\n"

        # Send the score vector to iB4e
        iB4e.stdin.write(result)

    build_sage_polytope_file(points, sagefile)
    
def build_sage_polytope_file(points, sagefile):
    pointstring = ", ".join(str(point) for point in points)
    file = open(sagefile, mode='w+')
    file.write("points = [" + pointstring + "]\n")
    file.write("poly = Polyhedron(points, base_ring=QQ)\n")
    file.write("slice = Polyhedron(eqns=[(1,0,0,0,1)])\n")
    file.write("fan = NormalFan(poly)\n")
    file.write("regions = [slice.intersection(cone.polyhedron()) for cone in fan.cones()[-1]]")
    file.close()

# Voodoo to make Python run the program
if __name__ == "__main__":
    main(sys.argv[1:])
