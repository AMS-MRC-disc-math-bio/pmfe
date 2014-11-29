// Copyright (c) 2014 Andrew Gainer-Dewar.

#ifndef IB4E_H
#define IB4E_H

#include <CGAL/Cartesian_d.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/Regular_complex_d.h>

#include <CGAL/Kernel_d/Point_d.h>
#include <CGAL/Kernel_d/Vector_d.h>
#include <CGAL/Kernel_d/Hyperplane_d.h>
#include <CGAL/Linear_algebraCd.h>

#include <CGAL/Origin.h>
#include <CGAL/Gmpq.h>

#include <string>

// Here come the typedefs!
typedef CGAL::Gmpq Q;

typedef CGAL::Cartesian_d<Q> K;
typedef CGAL::Point_d<K> QPoint;
typedef CGAL::Vector_d<K> QVector;
typedef CGAL::Hyperplane_d<K> Hyperplane;
typedef CGAL::Linear_algebraCd<Q> LinearAlgebra;
typedef LinearAlgebra::Matrix LMatrix;
typedef LinearAlgebra::Vector LVector;

typedef CGAL::Convex_hull_d<K> ConvexHull;
typedef ConvexHull::Facet_handle Facet;
typedef ConvexHull::Facet_iterator Facet_iterator;
typedef ConvexHull::Vertex_handle Vertex;
typedef ConvexHull::Vertex_iterator Vertex_iterator;

int iB4e_main(std::string seq_file_path, std::string param_dir = "Turner99");

#endif
