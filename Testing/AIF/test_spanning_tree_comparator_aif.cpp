// Copyright (c) 2012-2022 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of
// the License, or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

#include "FEVV/DataStructures/DataStructures_aif.h"

#include "FEVV/Wrappings/Geometry_traits_aif.h"
#include "FEVV/Wrappings/properties_aif.h"

#include "FEVV/Tools/Comparator/AIF/SpanningTreeVertexEdgeComparator.hpp" 

#include "FEVV/DataStructures/AIF/AIFMeshReader.hpp"
#include "FEVV/Tools/IO/FileUtilities.hpp"

#include <string>

using namespace FEVV;
using namespace FEVV::Comparator;

void
testSpanningTreeComparatorSurfaceMesh(const std::string& filename)
{
  using PtrMeshT = FEVV::DataStructures::AIF::AIFMesh::ptr_mesh;
  /**********************************************************************************************************/
  PtrMeshT pm;
  FEVV::DataStructures::AIF::AIFMeshReader in;
  if(!(pm = in.read(filename)))
  {
    std::cout << "Unable to read file " << filename << std::endl;
    return;
  }
  auto pos_pm = get(boost::vertex_point, *pm);
  /**********************************************************************************************************/
  auto st = get_spanning_tree_comparator(*pm, pos_pm, false);
}

int
main(int narg, char **argv)
{
  if(narg < 2)
  {
    std::cout << "Usage: " << argv[0]
              << " filename; filename being an off file." << std::endl;
    exit(EXIT_FAILURE);
  }

  testSpanningTreeComparatorSurfaceMesh(argv[1]);
  return 0;
}
