from sage.all import *
from sage.geometry.polyhedron.parent import Polyhedra
from sage.geometry.polyhedron.backend_ppl import Polyhedron_QQ_ppl
from sage.rings.rational_field import QQ
from sage.geometry.fan import Fan, Cone_of_fan
from collections import namedtuple
from pickle import UnpicklingError
import os.path

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
        """Return the normal fan of `self`."""
        return self._normal_fan

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

        points = []

        # Read the point data from the specified file
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

        filebase = os.path.splitext(polyfile)[0]

        # Construct the polytope
        try:
            # If the polytope is available in a pickle, we should use that
            thepoly = load(filebase) # Attempt to load the polytope from a pickle file

            # Verify that the vertices are the same as our points
            picklevertices = frozenset(map(tuple, thepoly.vertices()))
            filevertices = frozenset(tuple(point.vector) for point in points)
            if not picklevertices == filevertices:
                raise AssertionError

        except (IOError, UnpicklingError, AssertionError):
            # In any of these exception cases, the load failed, so we generate the polytope from the points
            thepoly = RNAPolytope(points)
            thepoly._build_normal_fan()
            thepoly.dump(filebase, -1)

        return thepoly
