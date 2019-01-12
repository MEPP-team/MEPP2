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
#include <iostream>
#include <fstream>
#include <OpenMesh/Core/IO/MeshIO.hh>

#define CGAL_USE_OM_POINTS
#include "FEVV/Wrappings/Geometry_traits_openmesh.h"

#include "FEVV/Filters/Generic/scaling.hpp"

// Boost properties adaptor
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/boost/graph/properties_PolyMesh_ArrayKernelT.h>

using namespace FEVV;
using namespace FEVV::Filters;

//------------------------------------------------------------------------------
void
test_calculate_scaling_open_mesh(char *filename)
{
  typedef OpenMesh::PolyMesh_ArrayKernelT<> Mesh;
  typedef FEVV::Geometry_traits< Mesh > Geometry;
  typedef Geometry::Scalar Scalar;

  std::ifstream in(filename);
  Mesh m;
  if(!OpenMesh::IO::read_mesh(m, filename))
  {
    std::cout << "failed";
    return;
  }

  auto pos_pm = get(boost::vertex_point, m);
  calculate_scaling(m, pos_pm, Scalar(2.0), Scalar(3.0), Scalar(4.0));

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

  test_calculate_scaling_open_mesh(argv[1]);
  return 0;
}
