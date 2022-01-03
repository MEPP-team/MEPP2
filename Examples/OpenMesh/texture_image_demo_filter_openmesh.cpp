// Copyright (c) 2012-2021 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published 
// by the Free Software Foundation; either version 3 of the License, 
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.



#include "Examples/Generic/TextureImageDemo/texture_image_demo_main.hpp"

// The following include defines the FEVV::MeshOpenMesh type
#include "FEVV/DataStructures/DataStructures_openmesh.h"

#include "FEVV/Wrappings/Geometry_traits_openmesh.h"
#include "FEVV/Wrappings/properties_openmesh.h"


// Main: load a mesh, apply the filter, write the mesh
int
main(int argc, const char **argv)
{
  return texture_image_demo_main< FEVV::MeshOpenMesh >(argc, argv);
}
