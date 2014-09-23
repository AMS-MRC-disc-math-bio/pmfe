from sage.all import *

def affine_hull(points, dim = None):
    from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace

    if dim == None:
        dim = len(points[0])

    ambient_space = VectorSpace(QQ, dim)
    vectors = [ambient_space(point) for point in points]

    shifted_vectors = [vector - vectors[0] for vector in vectors]
    shifted_linear_hull = ambient_space.subspace(shifted_vectors)

    the_affine_hull = AffineSubspace(vectors[0], shifted_linear_hull)
    return the_affine_hull

def facet_normal_vector(points, dim):
    ambient_space = VectorSpace(QQ, dim)
    vectors = [ambient_space(point) for point in points]
    center = Polyhedron(vertices = vectors).center()

    shifted_vectors = [vector - vectors[0] for vector in vectors]
    shifted_linear_hull = ambient_space.subspace(shifted_vectors)
    normal_space = shifted_linear_hull.complement()

    return normal_space.basis()[0]

def build_polytope(find_vector_oracle, dim):
    from itertools import chain
    ### INITIALIZE
    failed_target_vectors = []

    # Huggins' algorithm searches for the initial vectors using normals, but checking the basis set is simpler and should still catch everything
    ambient_space = VectorSpace(QQ, dim)
    initial_vector_pairs = ((find_vector_oracle(basis_element), find_vector_oracle(-basis_element)) for basis_element in ambient_space.basis())
    initial_vectors = (vector for pair in initial_vector_pairs for vector in pair)
    tentative_polytope = Polyhedron(vertices = initial_vectors)

    ### MAIN LOOP
    confirmed_facets = []
    while len(confirmed_facets) < tentative_polytope.n_facets():
        for index in range(tentative_polytope.n_facets()):
            facet = tentative_polytope.faces(dim-1)[index]
            vertex_vector_set = set(tuple(coord for coord in vert.vector()) for vert in facet.vertices())
            if vertex_vector_set in confirmed_facets:
                continue

            # TODO: Make sure this is finding the right normals
            innernormal = NormalFan(tentative_polytope).rays()[index]
            outernormal = -innernormal
            
            v = find_vector_oracle(outernormal)
            if v in facet.polyhedron():
                confirmed_facets.append(vertex_vector_set)
                continue
            else:
                vertexlist = list(tentative_polytope.vertices()) + [v]
                tentative_polytope = Polyhedron(vertices = vertexlist)
                break

    return tentative_polytope
