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
#include <CGAL/Cartesian.h>

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>

#include "Testing/Concepts/GeometryConceptCheck.h"
#include "FEVV/Wrappings/Geometry_traits_cgal_linear_cell_complex.h"

struct Myitem
{
  template< class Refs >
  struct Dart_wrapper
  {
    typedef CGAL::Dart< 2, Refs > Dart;
    typedef CGAL::Cell_attribute_with_point< Refs > Vertex_attribute;
    typedef CGAL::Cell_attribute< Refs > Face_attribute;
    typedef CGAL::cpp11::tuple< Vertex_attribute, void, Face_attribute >
        Attributes;
  };
};

int
main(int narg, char **argv)
{
  typedef CGAL::Cartesian< double > Kernel;
  typedef Kernel::Point_3 Point;
  typedef Kernel::Vector_3 Vector;
  typedef CGAL::Linear_cell_complex_traits< 3, Kernel > MyTraits;
  typedef CGAL::
      Linear_cell_complex_for_combinatorial_map< 2, 3, MyTraits, Myitem >
          Mesh;

  BOOST_CONCEPT_ASSERT((FEVV::Point_3Concept< Mesh >));

  typedef FEVV::Geometry_traits< Mesh > Geometry;
  Mesh mesh;
  Geometry g(mesh);
  // We don't have to call GeometryTraits::verify_assumption() because
  // we known that a CGAL::Vector_3 hardwires an ambiant space
  // dimension of 3:
  // FEVV::GeometryTraits::verify_assumption( g );
  FEVV::GeometryTraits::verify_operator_results(g);

  return 0;
}
