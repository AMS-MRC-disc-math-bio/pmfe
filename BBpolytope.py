from sage.all import *

polytope_facet_set = lambda polytope: set(frozenset(tuple(coord for coord in vert.vector()) for vert in facet.vertices()) for facet in polytope.faces(polytope.dim()-1))

def build_polytope(find_vector_oracle, dim):
    from itertools import chain
    ### INITIALIZE
    # Huggins' algorithm searches for the initial vectors using normals, which is algorithmically complicated in Sage.
    # Instead, we just check the 8 positive and negative basis vectors, which will still find us something full-dimensional.
    ambient_space = VectorSpace(QQ, dim)
    initial_vectors = [find_vector_oracle(basis_element) for basis_element in ambient_space.basis()] + [find_vector_oracle(-basis_element) for basis_element in ambient_space.basis()]
    tentative_polytope = Polyhedron(vertices = initial_vectors)

    ### MAIN LOOP
    confirmed_facets = set()

    # Keep looping until we have confirmed every face of the polytope.
    while len(confirmed_facets) < tentative_polytope.n_facets():
        # Test every facet of the polytope.
        for facet in tentative_polytope.faces(dim-1):
            # Test whether the facet is already confirmed.
            vertex_vector_set = frozenset(tuple(coord for coord in vert.vector()) for vert in facet.vertices())
            if vertex_vector_set in confirmed_facets:
                # If so, continue to the next facet.
                continue 

            # If the facet is not already confirmed, we need to test it.
            # We use the geometric center of the facet as the test vector.
            facet_center = facet.as_polyhedron().center()

            # We also need to find the outer normal of the facet.
            # Unfortunately, Sage's NormalFan method doesn't give these rays in a useful order, so we have to do it ourselves.
            # This is a nasty hack, but it works.
            # The ambient_Hrepresentation()[0] call gives us a representation of the facet as an inequality corresponding to a half-plane, and A() gives us the coefficient vector.
            # This is constructed in such a way that the result is the inner normal to the facet.
            # (If the facet is valid, every point in it will minimize dot product with this vector.)
            innernormal = ambient_space(facet.ambient_Hrepresentation()[0].A())

            # Use the oracle to find a vector minimizing with respect to the inner normal.
            v = ambient_space(find_vector_oracle(innernormal))
            # The heart of the algorithm.
            # If the dot product of the normal and v is less than that of the normal and the facet's center,
            # the facet is invalid.
            if innernormal.dot_product(v) < innernormal.dot_product(facet_center) and v not in tentative_polytope:
                # If so, add the vector to the polytope
                vertexlist = list(tentative_polytope.vertices()) + [v]
                tentative_polytope = Polyhedron(vertices = vertexlist)
                 # and discard any spurious confirmed facets using set intersection.
                confirmed_facets = confirmed_facets & polytope_facet_set(tentative_polytope)
                # Since we changed tentative_polytope, we return to the outer loop and begin checking facets again.
                break
            else:
                # Otherwise, confirm the facet
                confirmed_facets.add(vertex_vector_set)
                # and continue through the facet list.
                continue

    # If we made it out of the while loop, it should mean that we have confirmed every facet.
    # We test to make sure that's true.
    assert polytope_facet_set(tentative_polytope) <= confirmed_facets, "Not all facets confirmed! This is a serious algorithm error. Please contact the author with information about the sequence file that you are using."

    # and then return the result.
    return tentative_polytope
