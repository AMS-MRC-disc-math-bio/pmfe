from sage.all import *

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
    confirmed_facets = set()
    while len(confirmed_facets) < tentative_polytope.n_facets():
        for facet in tentative_polytope.faces(dim-1):
            vertex_vector_set = frozenset(tuple(coord for coord in vert.vector()) for vert in facet.vertices())
            if vertex_vector_set in confirmed_facets:
                continue

            facet_center = facet.as_polyhedron().center()

            #Ugly hack, but the NormalFan doesn't come out in a useful order
            normal = ambient_space(facet.ambient_Hrepresentation()[0][1:])
            
            v = ambient_space(find_vector_oracle(normal))
            if normal.dot_product(v) < normal.dot_product(facet_center) and v not in tentative_polytope:
                vertexlist = list(tentative_polytope.vertices()) + [v]
                tentative_polytope = Polyhedron(vertices = vertexlist)
                break
            else:
                confirmed_facets.add(vertex_vector_set)
                continue

    return tentative_polytope
