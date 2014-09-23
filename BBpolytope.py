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

    shifted_vectors = [vector - vectors[0] for vector in vectors]
    shifted_linear_hull = ambient_space.subspace(shifted_vectors)
    normal_space = shifted_linear_hull.complement()

    return normal_space.basis()[0]

def build_polytope(find_vector_oracle, dim):
    ### INITIALIZE
    polytope_vertices = []
    failed_target_vectors = []

    ambient_space = VectorSpace(QQ, dim)
    vertices_affine_hull = []

    ### Determine the affine hull of the polytope by testing vectors
    for k in xrange(dim+1):
        if len(polytope_vertices) != 0:
            vertices_affine_hull = affine_hull(polytope_vertices, dim)
            union_space = ambient_space.subspace(vertices_affine_hull.linear_part().basis() + failed_target_vectors)

        else:
            vertices_affine_hull = []
            union_space = ambient_space.subspace(failed_target_vectors)

        target_vector = union_space.complement().basis()[0]

        a = find_vector_oracle(target_vector)
        if a not in vertices_affine_hull:
            polytope_vertices.append(a)
            continue

        b = find_vector_oracle(-target_vector)
        if b not in vertices_affine_hull:
            polytope_vertices.append(b)
            continue

        failed_target_vectors.append(c)

    ### MAIN LOOP
    tentative_polytope = Polyhedron(vertices = polytope_vertices)

    confirmed_facets = []
    while len(confirmed_facets) < tentative_polytope.n_facets():
        for facet in tentative_polytope.faces(dim-1):
            vertex_vector_set = set(tuple(coord for coord in vert.vector()) for vert in facet.vertices())
            if vertex_vector_set in confirmed_facets:
                continue

            outernormal = facet_normal_vector(vertex_vector_set, dim)
            v = find_vector_oracle(outernormal)
            if v in facet.polyhedron():
                confirmed_facets.append(vertex_vector_set)
                continue
            else:
                polytope_vertices.append(v)
                tentative_polytope = Polyhedron(vertices = polytope_vertices)
                break

    return tentative_polytope
