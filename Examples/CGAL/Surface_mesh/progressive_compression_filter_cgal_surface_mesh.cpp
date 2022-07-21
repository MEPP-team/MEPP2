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

#include "FEVV/Filters/CGAL/Progressive_Compression/Helpers/initializer_compression.h"


// Main: load a mesh, apply the filter, write the mesh

//Input:
// path to a mesh
// from here, the parameters are optional, and have default values
// mode (0:compressing the mesh with defined parameters 
//       1: compressing the given mesh with every method combination. 
//       2: given a folder, will compute the distorsion between the meshes in the folder and the given mesh)
// path to the compressed mesh (binary file)
// folder to store the measure files 
// predictor as an int (see parameters.h for actual values)
// metric as an int
// kept position as an int
// number of simplification batches
// minimum number of vertices
// batch stopping condition as an int (0=ALL_EDGES or 1=REACH_THRESHOLD)
// quantization bits


//Output: simplified mesh (as progressive_compressionFilteroutput.obj )
//Compressed mesh as binary file(name given by the user or [predictor][keptposition][metric][quantization].bin
int
main(int argc, const char **argv)
{
  return progressive_compression_main< FEVV::MeshSurface >(argc, argv);
}
