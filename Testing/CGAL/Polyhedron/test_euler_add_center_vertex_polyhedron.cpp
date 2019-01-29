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
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Cartesian.h>
#include <CGAL/boost/graph/Euler_operations.h>

#include "FEVV/Wrappings/Geometry_traits_cgal_polyhedron_3.h"
#include "Testing/Utils/utils_are_meshes_identical.hpp"

#include <fstream>
#include <iostream>

using namespace FEVV;

//------------------------------------------------------------------------------

const char *outputFileName = "add_center_vertex_polyhedron_output.off";

//------------------------------------------------------------------------------

void
test_euler_add_center_vertex_polyhedron(std::string filename)
{
  typedef CGAL::Cartesian< double > Kernel;
  typedef CGAL::Polyhedron_3< Kernel, CGAL::Polyhedron_items_with_id_3 > Mesh;

  // load mesh

  std::ifstream in(filename);
  Mesh m;
  in >> m;

  // divise the first face, change topology only

  typedef boost::graph_traits< Mesh > GraphTraits;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
  typedef typename GraphTraits::face_descriptor face_descriptor;

  face_descriptor f = *(faces(m).begin());
  halfedge_descriptor h = halfedge(f, m);
  if(h == GraphTraits::null_halfedge())
  {
    std::cout << "Failed to retrieve the first face halfedge. Exiting."
              << std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout << "Adding center vertex..." << std::endl;
  halfedge_descriptor hnew = CGAL::Euler::add_center_vertex(h, m);

  // fix geometry

  typedef GraphTraits::vertex_descriptor vertex_descriptor;
  vertex_descriptor vnew = target(hnew, m);

  typedef FEVV::Geometry_traits< Mesh > Geometry;
  typedef Geometry::Point Point;
  put(get(boost::vertex_point, m), vnew, Point(4.0f, 3.0f, 0.0f));

  // Write result to file

  std::ofstream out(outputFileName);
  out << m;
}

//------------------------------------------------------------------------------

int
main(int narg, char **argv)
{
  if(narg < 2)
  {
    std::cout << "Usage: " << argv[0] << "  mesh_file  [reference_result_file]"
              << std::endl;
    std::cout << "Example: " << argv[0] << "  airplane.off" << std::endl;
    exit(EXIT_FAILURE);
  }

  // run test
  test_euler_add_center_vertex_polyhedron(argv[1]);

  // compare with reference result
  if(narg == 3)
  {
    if(!are_meshes_equal(outputFileName, argv[2], false))
      return 1; // test failed
  }

  return 0;
}
