// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published 
// by the Free Software Foundation; either version 3 of the License, 
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#include <iostream>
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include "FEVV/Wrappings/Geometry_traits_cgal_polyhedron_3.h"

int
main(int narg, char **argv)
{
  if(narg > 1)
  {
    std::cout << "Usage: " << argv[0] << " (no arguments)." << std::endl;
    exit(EXIT_FAILURE);
  }

  typedef CGAL::Cartesian< double > Kernel;
  typedef CGAL::Polyhedron_3< Kernel, CGAL::Polyhedron_items_with_id_3 > Mesh;

  typedef FEVV::Geometry_traits< Mesh > Geometry;
  typedef Geometry::Point Point;
  typedef Geometry::Vector Vector;

  Mesh m;
  Geometry g(m);
  Point p1(0.0, 0.0, 0.0);
  Point p2(1.0, 0.0, 0.0);
  Point p3(0.0, 0.1, 0.0);

  Vector n = g.unit_normal(p1, p2, p3);

  if(n == Vector(0.0, 0.0, 1.0))
  {
    std::cout << "OK." << std::endl;
    return 0;
  }
  std::cout << "Result is " << n << " when (0, 0, 1) was expected."
            << std::endl;
  return 1;
}
