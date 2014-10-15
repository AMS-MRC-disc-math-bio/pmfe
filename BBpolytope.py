from sage.all import *

facet_vertex_set = lambda facet: frozenset(tuple(coord for coord in vert.vector()) for vert in facet.vertices())
polytope_facet_set = lambda polytope: set(facet_vertex_set(facet) for facet in polytope.faces(polytope.dim()-1))

def build_polytope(find_vector_oracle, dim, maximize=True):
    from itertools import chain
    ### INITIALIZE
    # Huggins' algorithm searches for the initial vectors using normals, which is algorithmically complicated in Sage.
    # Instead, we just check the 8 positive and negative basis vectors, which will still find us something full-dimensional.
    ambient_space = VectorSpace(QQ, dim)
    initial_vectors = [find_vector_oracle(basis_element) for basis_element in ambient_space.basis()] + [find_vector_oracle(-basis_element) for basis_element in ambient_space.basis()]
    tentative_polytope = Polyhedron(vertices = initial_vectors)

    # Sanity check: make sure that the polytope is full-dimensional at this point.
    assert tentative_polytope.dimension() == dim, "Initial polytope is not full-dimensional!"

    ### MAIN LOOP
    confirmed_facets = set()

    # Keep looping until we have confirmed every face of the polytope.
    while not polytope_facet_set(tentative_polytope).issubset(confirmed_facets):
        # Test every facet of the polytope.
        for facet in tentative_polytope.faces(dim-1):
            # Test whether the facet is already confirmed.
            vertex_vector_set = facet_vertex_set(facet)
            if vertex_vector_set in confirmed_facets:
                # If so, continue to the next facet.
                continue 

            # If the facet is not already confirmed, we need to test it.
            # To do this, we need to find the inner normal of the facet.
            # Unfortunately, Sage's NormalFan method doesn't give these rays in a useful order, so we have to do it ourselves.
            # This is a nasty hack, but it works.
            # The ambient_Hrepresentation()[0] call gives us a representation of the facet as an inequality corresponding to a half-plane, and A() gives us the coefficient vector.
            # (If the facet is valid, every point in it will minimize dot product with this vector.)
            innernormal = ambient_space(facet.ambient_Hrepresentation()[0].A())
            testnormal = -innernormal if maximize else innernormal
            facetscore = -facet.ambient_Hrepresentation()[0].b()

            # Use the oracle to find an optimal vector
            v = ambient_space(find_vector_oracle(testnormal))
            # To test the facet, we compare the dot product of the testnormal with v.
            # Which test we do depends on whether the underlying oracle is a maximizer or minimizer.
            if maximize and testnormal.dot_product(v) > facetscore or not maximize and testnormal.dot_product(v) < facetscore:
                # If so, add the vector to the polytope
                vertexlist = list(tentative_polytope.vertices()) + [v]
                new_tentative_polytope = Polyhedron(vertices = vertexlist)
                assert new_tentative_polytope != tentative_polytope, "Facet confirmation error!"
                tentative_polytope = new_tentative_polytope
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
