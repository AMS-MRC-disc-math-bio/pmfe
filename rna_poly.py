from sage.all import *
from sage.geometry.polyhedron.parent import Polyhedra
from sage.geometry.polyhedron.backend_ppl import Polyhedron_QQ_ppl
from sage.rings.rational_field import QQ
from sage.geometry.fan import Fan, Cone_of_fan
from collections import namedtuple
from pickle import UnpicklingError
import os.path

pickle_version = 2

class RNAPolytope(Polyhedron_QQ_ppl):
    def __init__(self, points):
        """
        Construct an RNAPolytope from a collection of points.

        Keyword arguments:
        points -- An Iterable of points, each of which has (at least) members `structure` (holding arbitrary data about the point) and `vector` (giving the coordinates of the point in a format which can be handled by `Polytope`)
        """
        self._normal_fan = None
        self._vertex_from_cone = None

        parent = Polyhedra(QQ, len(points[0].vector))

        vertex_list = [point.vector for point in points]
        super(RNAPolytope, self).__init__(parent, [vertex_list, None, None], None)

        structure_dict = {tuple(map(QQ, point.vector)): point.structure for point in points}
        self._structures = structure_dict

    def _build_normal_fan(self):
        if self._normal_fan is None:
            cones = [[ieq.index() for ieq in vertex.incident()] for vertex in self.vertices()]
            rays = [ ieq.A() for ieq in self.inequalities() ]

            fan = Fan(cones, rays, check=False, is_complete = True)
            self._normal_fan = fan

        if self._vertex_from_cone is None:
            conedict = {}
            for vertex in self.vertices():
                cone = self.normal_fan().cone_containing([ieq.A() for ieq in vertex.incident()])
                conedict[cone] = vertex
            self._vertex_from_cone = conedict

    def normal_fan(self):
        """
        Return the normal fan of `self`.
        """
        if not hasattr(self, "_normal_fan"):
            self._build_normal_fan()

        return self._normal_fan

    def _build_d1_slices(self):
        subspace = Polyhedron(eqns=[(-1, 0, 0, 0, 1)]) # Hyperplane d = 1
        slices = set()

        # Construct the slices
        for cone in self.normal_fan().cones(codim=0):
            coneslice = cone.polyhedron().intersection(subspace)

            # If the slice does not have volume, we're done
            if coneslice.dimension() < 3:
                continue

            # Otherwise, project to Q3
            vecprojector = lambda vec: vec[0:3] # Drop the d coordinate
            newverts = lambda polyslice: [vecprojector(vert) for vert in polyslice.vertices()] # Projects the defining vertices of a slice into Q3
            newrays = lambda polyslice: [vecprojector(ray) for ray in polyslice.rays()] # Projects the defining rays of a slice into Q3
            polyprojector = lambda polyslice: Polyhedron(vertices = newverts(polyslice), rays = newrays(polyslice)) # Projects a polytope from Q4 to Q3 using the helper functions we just defined

            # We can then apply polyprojector to our slice to move it down to Q3
            projslice = polyprojector(coneslice)

            # For future use, we'll store the original cone and the slice in Q4 as a member of this object
            projslice.original_cone = cone
            projslice.original_slice = coneslice

            # As well as the original vertex and structure associated to those cones
            projslice.original_vertex = self.vertex_from_cone(cone)
            projslice.original_structure = self[projslice.original_vertex]

            # And, finally, we add this slice to our list
            slices.add(projslice)

        # Once we've finished the loop, all the slices have been processed, so we store the results
        self._d1_slices = slices

    def d1_slices(self):
        """
        Return the d=1 slices of self as Polyhedra in Q3
        """

        if not hasattr(self, "_d1_slices"):
            self._build_d1_slices()

        return self._d1_slices

    def vertex_from_cone(self, cone):
        assert cone in self.normal_fan().cones(self.dim())
        return self._vertex_from_cone[cone]

    def __getitem__(self, key):
        try:
            processed_key = tuple(map(QQ, key)) # Try to cast the input to a tuple of rationals
            return self._structures[processed_key] # Cast the input to a tuple of rationals
        except TypeError:
            return None # Bad key!

    @classmethod
    def construct_from_file(cls, polyfile):
        """
        Construct an RNAPolytope from a file in polyfile format.

        Keyword arguments:
        polyfile -- A file representing the polytope structure
        """

        filebase = os.path.splitext(polyfile)[0]

        # Construct the polytope
        try:
            # If the polytope is available in a pickle, we should use that
            thepoly = load(filebase) # Attempt to load the polytope from a pickle file

            # If the pickled polytope is obsolete, however, we need to rebuild it
            if cls.poly_is_obsolete(thepoly):
                raise UnpicklingError

        except (IOError, UnpicklingError, AssertionError):
            # In any of these exception cases, the load failed, so we generate the polytope from the points

            # Read the point data from the specified file
            points = []
            with open(polyfile) as f:
                for line in f:
                    newline = line.split('#')[0] # Only use content before the comment symbol
                    if newline.strip() != '':
                        elts = newline.split()
                        structure = elts[1]
                        multiloops = QQ(elts[2])
                        unpaired = QQ(elts[3])
                        branches = QQ(elts[4])
                        w = QQ(elts[5])
                        energy = QQ(elts[6])
                        coords = (multiloops, unpaired, branches, w)

                        newpoint = namedtuple('RNApoint', ['structure', 'vector', 'energy'])
                        newpoint.vector = coords
                        newpoint.structure = structure
                        newpoint.energy = energy

                        points.append(newpoint)

            thepoly = RNAPolytope(points)
            thepoly._pickle_version = pickle_version
            thepoly._build_normal_fan()
            thepoly._build_d1_slices()
            thepoly.dump(filebase, -1)

        return thepoly

    @classmethod
    def poly_is_obsolete(cls, poly):
        """
        Test whether the polytope object is obsolete, in which case it should be upgraded
        """

        if (not hasattr(poly, "_pickle_version") or poly._pickle_version < pickle_version):
            return true
        else:
            return false
