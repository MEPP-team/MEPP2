// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of
// the License, or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#include <iostream>
#include <fstream>

#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>
#include "FEVV/Wrappings/Geometry_traits_cgal_linear_cell_complex.h"

#include <CGAL/Linear_cell_complex_incremental_builder.h>
#include <CGAL/Combinatorial_map_constructors.h>

#include <CGAL/Cartesian.h>

typedef CGAL::Cartesian< double > Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Linear_cell_complex_traits< 3, Kernel > MyTraits;

struct Myitem
{
  template< class Refs >
  struct Dart_wrapper
  {
    typedef CGAL::Tag_true Darts_with_id;
    typedef CGAL::Cell_attribute_with_point_and_id< Refs > Vertex_attribute;
    typedef CGAL::Cell_attribute_with_id< Refs > Face_attribute;
    typedef CGAL::cpp11::tuple< Vertex_attribute, void, Face_attribute >
        Attributes;
  };
};

typedef CGAL::
    Linear_cell_complex_for_combinatorial_map< 2, 3, MyTraits, Myitem >
        LCC;

//------------------------------------------------------------------------------
int
main(int narg, char **argv)
{
  LCC lcc;
  typedef LCC::Dart_handle Dart_handle;

  try
  {
#if 1
    Dart_handle s1 = lcc.make_segment(
        Point(1, 0, 10), Point(0, 0, 10)); // add a dangling edge (non-manifold)

    Dart_handle t1 =
        lcc.make_triangle(Point(1, 0, 1), Point(0, 0, 0), Point(1, 1, 0));
    Dart_handle t2 =
        lcc.make_triangle(Point(1, 1, 1), Point(0, 0, 0), Point(1, 1, 0));
    Dart_handle t3 =
        lcc.make_triangle(Point(1, -1, 1), Point(0, 0, 0), Point(1, 1, 0));
#else
    Dart_handle s1 = lcc.make_edge();
    Dart_handle t1 = lcc.make_combinatorial_polygon(3);
    Dart_handle t2 = lcc.make_combinatorial_polygon(3);
    Dart_handle t3 = lcc.make_combinatorial_polygon(3);
#endif
    // std::cout << "Before sewing: triangle 1, p1: " <<
    // lcc.attribute<0>(t1)->point() << std::endl; std::cout << "Before sewing:
    // triangle 1, p2: " << lcc.attribute<0>(lcc.beta(t1, 1))->point() <<
    // std::endl; std::cout << "Before sewing: triangle 1, p3: " <<
    // lcc.attribute<0>(lcc.beta(t1, 0))->point() << std::endl;
    lcc.sew< 2 >(lcc.beta(t1, 1), lcc.beta(t2, 1));
    // std::cout << "After sewing: triangle 1, p1: " <<
    // lcc.attribute<0>(t1)->point() << std::endl; std::cout << "After sewing:
    // triangle 1, p2: " << lcc.attribute<0>(lcc.beta(t1, 1))->point() <<
    // std::endl; std::cout << "After sewing: triangle 1, p3: " <<
    // lcc.attribute<0>(lcc.beta(t1, 0))->point() << std::endl;

    // lcc.display_characteristics(std::cout);
    // std::cout << std::endl;

    lcc.sew< 2 >(
        lcc.beta(t2, 1),
        lcc.beta(
            t3,
            1)); // try to sew 3 triangles around the same edge => exception,
                 // because LCC cannot handle non-manifold topologies
  }
  catch(...)
  {
    std::cerr << "test_not_2_manifold_linear_cell_complex: LCC and "
                 "combinatorial maps cannot handle non-manifold topologies"
              << std::endl;
  }

  return EXIT_SUCCESS;
}
