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
#if(_MSC_VER >= 1400)
#ifndef _SCL_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif
#endif

#include <OpenMesh/Core/IO/MeshIO.hh>
// Graph traits adaptors
#define CGAL_USE_OM_POINTS
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>

#include <fstream>
#include <iostream>

#include "Testing/Utils/utils_retrieve_halfedge.h"
#include "Testing/Utils/utils_are_meshes_identical.hpp"
#include "FEVV/Operators/Generic/Manifold/collapse_edge_euler.hpp"
// DBG #include "FEVV/Filters/Generic/print_points.hpp"

using namespace FEVV;
// using namespace FEVV::Filters;

//------------------------------------------------------------------------------

void
test_collapse_edge_open_mesh(std::string filename,
                             int source_index,
                             int target_index,
                             const std::string &output_file_name)
{
  typedef OpenMesh::PolyMesh_ArrayKernelT< /*MyTraits*/ > Mesh;

  // Load a mesh
  std::ifstream in(filename);
  Mesh m;
  if(!OpenMesh::IO::read_mesh(m, filename))
  {
    std::cout << "failed to open file " << filename << "\n";
    exit(EXIT_FAILURE);
  }

  // Collapse the first edge
  typedef boost::graph_traits< Mesh > GraphTraits;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;

  halfedge_descriptor h = retrieve_halfedge(m, source_index, target_index);

  if(h == GraphTraits::null_halfedge())
  {
    std::cout << "Failed to retrieve edge from " << source_index << " to "
              << target_index << "." << std::endl;
    std::cout << "Exiting";
    exit(EXIT_FAILURE);
  }

  std::cout << "Collapsing edge " << source_index << " to " << target_index
            << "." << std::endl;
  // DBG FEVV::Filters::print_points(m, get(CGAL::vertex_point, m));
  FEVV::Operators::collapse_edge_keep_target_euler(m, h);
  // DBG FEVV::Filters::print_points(m, get(CGAL::vertex_point, m));

  // Write result to file
  m.garbage_collection();
  // Required for the geometry to be cleaned before
  // writing to file ; maybe this is a bug in OpenMesh
  OpenMesh::IO::write_mesh(m, output_file_name);
}

//------------------------------------------------------------------------------

int
main(int narg, char **argv)
{
  if(narg < 3 || narg > 4)
  /// \todo This should be the same as test_collapse_edge_polyhedron.cpp
  {
    std::cout << "Usage: a.out filename; filename being an off file."
              << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string to_treat = argv[1];

  int deletion_case = std::stoi(std::string(argv[2]));
  std::string output_file_name = std::string(argv[0]) + argv[2] + ".off";

  if(deletion_case == 0)
  {
    test_collapse_edge_open_mesh(
        to_treat, 0, 1, output_file_name); // Whatever comes first
  }
  else if(deletion_case == 1)
  {
    test_collapse_edge_open_mesh(to_treat, 3, 4, output_file_name);
  }
  else if(deletion_case == 2)
  {
    test_collapse_edge_open_mesh(to_treat, 4, 3, output_file_name);
  }

  if(narg == 4)
  {
    if(!are_meshes_equal(output_file_name, argv[3], false))
      return 1; // test failed
  }

  return 0; // test succeeded
}
