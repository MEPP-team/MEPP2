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

#include "FEVV/Wrappings/Geometry_traits_aif.h"
#include "FEVV/Wrappings/properties_aif.h"

/// \todo Interestingly enough the following inclusion must be placed _after_
/// the AIF specific include files for the file to compile.
/// This behaviour is different from the ones of the equivalent examples
/// using other data structures like OpenMesh or the ones of CGAL.
/// Inquire on this (peculiar?) behavioral difference of AIF.
#include "Examples/Generic/Helloworld/helloworld_main.hpp"

int
main(int argc, const char **argv)
{
  return helloworld_main< FEVV::MeshAIF >(argc, argv);
}
