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
#include <fstream>

// Surface mesh
#include <CGAL/Surface_mesh.h>
#include <CGAL/Cartesian.h>

// Graph traits adaptors
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include "FEVV/Wrappings/Geometry_traits_cgal_surface_mesh.h"

#include "FEVV/Filters/Generic/calculate_face_normals.hpp"
#include "FEVV/Filters/Generic/scaling.hpp"
#include "FEVV/Operators/Generic/Manifold/collapse_edge.hpp"
#include "FEVV/Filters/Generic/print_points.hpp"
#include "FEVV/Filters/Generic/print_face_normals.hpp"

using namespace FEVV;
using namespace FEVV::Filters;

//------------------------------------------------------------------------------
void
test_surface_mesh(char *filename)
{
  typedef CGAL::Cartesian< double > Kernel;
  typedef Kernel::Point_3 Point;
  typedef Kernel::Vector_3 Vector;
  typedef CGAL::Surface_mesh< Point > Mesh;

  // 1) Load a mesh
  std::ifstream in(filename);
  Mesh m;
  try
  {
    in >> m;
  }
  catch(const std::length_error &le)
  {
    std::cerr << "[SurfaceMesh] Exception caught while reading input file "
              << filename << ": " << le.what() << std::endl;
    BOOST_ASSERT_MSG(
        false, "[SurfaceMesh] Exception caught while reading input file.");
  }

  auto pos_pm = get(CGAL::vertex_point, m);

  // 2) Add face normal property
  Mesh::Property_map< Mesh::Face_index, Vector > normals;
  normals = m.add_property_map< Mesh::Face_index, Vector >("f:Normal").first;

  // 3) Compute normals
  calculate_face_normals(m // Graph
                         ,
                         pos_pm // map from vertex_descriptor to point
                         ,
                         normals // map from face_descriptor to Vector_3
  );

  // 4) Draw normals
  print_face_normals(m, normals);

  // 5) Scale the mesh
  calculate_scaling(m, pos_pm, 2.0, 3.0, 4.0);

  // 6) Display the points
  print_points(m, pos_pm);

  // 7) Collapse some edges
  /*std::cout << "Start collapsing edges ..." << std::endl;
  collapse_some_edges(M);
  std::ofstream collapse_some_edges_os("collapse_some_edges.off");
  collapse_some_edges_os << M;*/

  std::cout << "Done." << std::endl;
}

//------------------------------------------------------------------------------
int
main(int narg, char **argv)
{
  if(narg < 2)
  {
    std::cout
        << "Usage: test-surface-mesh filename; filename being an off file."
        << std::endl;
    exit(EXIT_FAILURE);
  }

  test_surface_mesh(argv[1]);
  return 0;
}
