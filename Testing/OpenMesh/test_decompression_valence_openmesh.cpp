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
#include "FEVV/DataStructures/DataStructures_openmesh.h"

#include "FEVV/Wrappings/Geometry_traits_openmesh.h"
#include "FEVV/Wrappings/properties_openmesh.h"

#include "Testing/Generic/Manifold/test_decompression_valence.inl"


int
main(int argc, const char **argv)
{
  return test_decompression_valence< FEVV::MeshOpenMesh >(argc, argv);
}
