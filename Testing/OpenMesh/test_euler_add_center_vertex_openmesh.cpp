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
#include "FEVV/Wrappings/Geometry_traits_openmesh.h" // Always in first position for MVS 4996 warnings

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include <CGAL/boost/graph/Euler_operations.h>

// OpenMesh adaptors
#define CGAL_USE_OM_POINTS
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/boost/graph/properties_PolyMesh_ArrayKernelT.h>

#include "Testing/Utils/utils_are_meshes_identical.hpp"

#include <fstream>
#include <iostream>

using namespace FEVV;

//------------------------------------------------------------------------------

const char *outputFileName = "add_center_vertex_openmesh_output.off";

//------------------------------------------------------------------------------

void
test_euler_add_center_vertex_open_mesh(std::string filename)
{
  typedef OpenMesh::PolyMesh_ArrayKernelT<> Mesh;

  // load mesh

  Mesh m;
  if(!OpenMesh::IO::read_mesh(m, filename))
  {
    std::cout << "failed to open file " << filename << "\n";
    exit(EXIT_FAILURE);
  }

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

  m.garbage_collection();
  // Required for the geometry to be cleaned before
  // writing to file ; maybe this is a bug in OpenMesh
  OpenMesh::IO::write_mesh(m, outputFileName);
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
  test_euler_add_center_vertex_open_mesh(argv[1]);

  // compare with reference result
  if(narg == 3)
  {
    if(!are_meshes_equal(outputFileName, argv[2], false))
      return 1; // test failed
  }

  return 0;
}
