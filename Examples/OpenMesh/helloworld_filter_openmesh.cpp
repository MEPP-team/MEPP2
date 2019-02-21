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
// The following include is generic: it works without any change for
// all mesh datastructures. This should be true for all code to be
// encountered in directories mentioning "Generic" within their pathname.
#include "Examples/Generic/Helloworld/helloworld_main.hpp"

// The following include defines the FEVV::MeshOpenMesh type that is
// used below when instantiating helloworld_main<>
#include "FEVV/DataStructures/DataStructures_openmesh.h"

// The filter itself (Examples/Generic/Helloworld/helloworld_filter.hpp),
// as included by helloworld_main.hpp, manipulates both the geometry
// of the mesh as well as some of the mesh properties. The following
// includes respectively allow the filter to use expressions like
// `FEVV::Geometry_traits<...>` and
// `FEVV::Face_pmap<HalfedgeGraph, float> my_property_map;`.
#include "FEVV/Wrappings/Geometry_traits_openmesh.h"
#include "FEVV/Wrappings/properties_openmesh.h"


// Main: load a mesh, apply the filter, write the mesh
int
main(int argc, const char **argv)
{
  return helloworld_main< FEVV::MeshOpenMesh >(argc, argv);
}
