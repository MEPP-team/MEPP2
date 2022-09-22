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

// The following include defines the FEVV::MeshSurface type
#include "FEVV/DataStructures/DataStructures_cgal_surface_mesh.h"

#include "FEVV/Wrappings/Geometry_traits_cgal_surface_mesh.h"
#include "FEVV/Wrappings/properties_surface_mesh.h"

#include "FEVV/Filters/CGAL/Progressive_Compression/Helpers/initializer_decompression.h"


// Main: load a compressed mesh file (binary file), apply the filter, write the reconstructed mesh

//Input:
// path to the compressed mesh (binary file)
// path to the output mesh filename: set by default to progressive_decompression_filter_output.obj
// number of maximum refinement batches to reconstruct: set by default to 10000

//Output: reconstructed/decompressed mesh
int
main(int argc, const char **argv)
{
  return progressive_decompression_main< FEVV::MeshSurface >(argc, argv);
}
