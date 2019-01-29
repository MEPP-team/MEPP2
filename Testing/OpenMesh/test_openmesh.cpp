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

#include <fstream>

#include "FEVV/Filters/Generic/calculate_face_normals.hpp"
#include "FEVV/Filters/Generic/print_face_normals.hpp"

// OpenMesh
#include "FEVV/DataStructures/DataStructures_openmesh.h"

using namespace FEVV;
using namespace FEVV::Filters;

//------------------------------------------------------------------------------
void
test_open_mesh(char *filename)
{
  typedef MeshOpenMesh Mesh;
  typedef Mesh::Normal Vector;

  // 1) Load a mesh
  std::ifstream in(filename);
  Mesh m;
  try
  {
    if(!OpenMesh::IO::read_mesh(m, filename))
    {
      std::cout << "failed";
      return;
    }
  }
  catch(const std::length_error &le)
  {
    std::cerr << "[OpenMesh] Exception catched while reading input file "
              << filename << ": " << le.what() << std::endl;
    BOOST_ASSERT_MSG(false,
                     "[OpenMesh] Exception catched while reading input file.");
  }

  // 2) Add face normal property
  typedef boost::property_map< Mesh, boost::face_index_t >::const_type
      Face_index_map;
  Face_index_map face_index_pm = get(boost::face_index, m);
  boost::vector_property_map< Vector, Face_index_map > normals(face_index_pm);

  // 3) Compute normals
  calculate_face_normals(
      m,                           // Graph
      get(boost::vertex_point, m), // map from vertex_descriptor to point
      normals                      // map from face_descriptor to Vector_3
  );

  // 4) Draw normals
  print_face_normals(m, normals);

  std::cout << "Done." << std::endl;
}

//------------------------------------------------------------------------------
int
main(int narg, char **argv)
{
  if(narg < 2)
  {
    std::cout << "Usage: " << argv[0]
              << " filename; filename being an off file." << std::endl;
    exit(EXIT_FAILURE);
  }

  test_open_mesh(argv[1]);
  return 0;
}
