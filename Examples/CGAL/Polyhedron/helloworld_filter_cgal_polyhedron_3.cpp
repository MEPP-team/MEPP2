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



#include "Examples/Generic/Helloworld/helloworld_main.hpp"

// The following include defines the FEVV::MeshPolyhedron type
#include "FEVV/DataStructures/DataStructures_cgal_polyhedron_3.h"

#include "FEVV/Wrappings/Geometry_traits_cgal_polyhedron_3.h"
#include "FEVV/Wrappings/properties_polyhedron_3.h"


// Main: load a mesh, apply the filter, write the mesh
int
main(int argc, const char **argv)
{
  return helloworld_main< FEVV::MeshPolyhedron >(argc, argv);
}
