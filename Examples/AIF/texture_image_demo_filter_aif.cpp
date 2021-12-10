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

// The following include defines the FEVV::MeshAIF type
#include "FEVV/DataStructures/AIF/AIFMesh.hpp"

#include "FEVV/Wrappings/Geometry_traits_aif.h"
#include "FEVV/Wrappings/properties_aif.h"

// Must be included AFTER AIF specific include files to avoid a
// conversion warning/error on Windows
// TODO inquire on this peculiar behavior of AIF
#include "Examples/Generic/TextureImageDemo/texture_image_demo_main.hpp"


// Main: load a mesh, apply the filter, write the mesh
int
main(int argc, const char **argv)
{
  return texture_image_demo_main< FEVV::MeshAIF >(argc, argv);
}
