#!/usr/bin/env python
import os, sys, argparse, subprocess
import GTrunner, GTsetMBparam, GTscorer

iB4e_path = "iB4e/iB4e-rna"

def main(argv):    
    # Set up parameters
    parser = argparse.ArgumentParser(description="Run GTfold with specified a, b, c, d parameters for iB4e")
    parser.add_argument("-s", "--sequence", nargs=1, help="Sequence to fold", required=True)

    args = vars(parser.parse_args())
    seqfile = args["sequence"][0]
    sagefile = os.path.splitext(seqfile)[0] + ".sage"

    turnerdir = "Turner99"
    outputdir = "output/data"
    scoredir = "output/scoring"

    points = []

    iB4e = subprocess.Popen([iB4e_path, '4'], # Our problem operates in 4 dimensions
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            shell=False)

    # Compute classical scores for later reference
    classical_struct = GTrunner.run_gt(turnerdir, "", seqfile)
    classical_scores = GTscorer.find_xyzw(turnerdir, scoredir, classical_struct)

    # Now turn iB4e loose to compute the geometry
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
        scores = GTscorer.find_xyzw(turnerdir, scoredir, structfile)
        points.append(scores)
        result = " ".join(map(str, scores)) + "\n"

        # Clean up the temporary structure file
        os.remove(structfile)

        # Send the score vector to iB4e
        iB4e.stdin.write(result)

    # Now build a Sage file encoding the desired data
    build_sage_polytope_file(classical_scores, points, sagefile)
    print "Wrote Sage polytope file " + sagefile
    
def build_sage_polytope_file(classical_scores, points, sagefile):
    pointstring = ", ".join(str(point) for point in points)
    file = open(sagefile, mode='w+')
    file.write("points = [" + pointstring + "]\n")
    file.write("classical_scores = " + str(classical_scores) + "\n")
    file.write("poly = Polyhedron(points, base_ring=QQ)\n")
    file.write("slice = Polyhedron(eqns=[(1,0,0,0,1)])\n")
    file.write("fan = NormalFan(poly)\n")
    file.write("regions4 = [slice.intersection(cone.polyhedron()) for cone in fan.cones()[-1]]\n")
    file.write("regions4 = filter(lambda region: not region.is_empty(), regions4)\n")
    file.write("vecprojector = lambda vec: vec[:-1]\n")
    file.write("polyprojector = lambda polyslice: Polyhedron(ieqs = [vecprojector(ieq) for ieq in polyslice.Hrepresentation()[1:]])\n")
    file.write("regions = [polyprojector(region) for region in regions4]\n")
    file.close()

# Voodoo to make Python run the program
if __name__ == "__main__":
    main(sys.argv[1:])
