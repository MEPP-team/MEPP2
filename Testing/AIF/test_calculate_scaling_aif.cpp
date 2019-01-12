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
#include "FEVV/DataStructures/AIF/AIFMesh.hpp"
#include "FEVV/Wrappings/Graph_traits_aif.h"
#include "FEVV/Wrappings/Geometry_traits_aif.h"

#include "FEVV/Wrappings/properties_aif.h"

using AIFKernel = FEVV::AIF_mesh_kernel_generator;
using MeshT = FEVV::DataStructures::AIF::AIFMesh;

//---------------- specific code above --------------------

//---------------- generic code below ---------------------

//#include "Testing/Utils/utils_retrieve_halfedge.h"
//#include "Testing/Utils/utils_are_meshes_identical.hpp"
#include "FEVV/Filters/Generic/scaling.hpp"

#include "FEVV/Filters/Generic/generic_reader.hpp"
// include "FEVV/Filters/Generic/generic_writer.hpp"

//---------------- generic code above ---------------------

//---------------- other code below -----------------------
#include <fstream>
#include <iostream>

void
test_calculate_scaling_aif(const std::string &input_file_path)
{
  typedef FEVV::Geometry_traits< MeshT > Geometry;
  typedef Geometry::Scalar Scalar;

  // read mesh from file
  MeshT m;
  FEVV::PMapsContainer pmaps_bag;
  FEVV::Filters::read_mesh(input_file_path, m, pmaps_bag);

  try
  {
    auto pos_pm = get(boost::vertex_point, m);
    FEVV::Filters::calculate_scaling(
        m, pos_pm, Scalar(2.0), Scalar(3.0), Scalar(4.0));
  }
  catch(...)
  {
    std::cout << "processing failed";
    return;
  }
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

  test_calculate_scaling_aif(argv[1]);
  return 0;
}
