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
#include <CGAL/Cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include "FEVV/Wrappings/Geometry_traits_cgal_surface_mesh.h"
#include "FEVV/Wrappings/Graph_traits_extension_cgal_surface_mesh.h"

#include "FEVV/Wrappings/properties_surface_mesh.h"
#include "FEVV/Filters/Generic/generic_reader.hpp"
#include "FEVV/Filters/Generic/generic_writer.hpp"
#include "FEVV/Filters/Generic/calculate_face_normals.hpp"

#include "FEVV/Tools/IO/FileUtilities.hpp"
#include "Testing/Utils/utils_are_meshes_identical.hpp"
#include "Testing/Utils/utils_identical_text_based_files.hpp"

#include "FEVV/Filters/CGAL/Surface_mesh/msdm2.h"
#include "FEVV/Filters/Generic/minmax_map.h"
#include "FEVV/Filters/Generic/color_mesh.h"

#include <iostream>
#include <string>


using CGALKernel = CGAL::Cartesian< double >;
using CGALPoint = CGALKernel::Point_3;
using MeshT = CGAL::Surface_mesh< CGALPoint >;


int
main(int argc, const char **argv)
{
  if(argc != 5)
  {
    std::cout
        << "Load 2 meshes from input files, compute their MSDM2, write it to "
           "an output file, then compare the output file with a reference file."
        << std::endl;
    std::cout << "Usage: " << argv[0]
              << "  input_original_mesh_file  input_deformed_mesh_file  "
                 "output_mesh_file  reference_mesh_file"
              << std::endl;
    std::cout << "Example: " << argv[0]
              << "  armadillo.obj armadillo_simplified.obj armadillo.out.obj "
                 "armadillo_output_reference.obj"
              << std::endl;
    return EXIT_FAILURE;
  }

  std::string input_original_file_path = argv[1];
  std::string input_deformed_file_path = argv[2];
  std::string output_file_path = argv[3];
  std::string reference_file_path = argv[4];

  //--------------------------------------------------

  int i = 0;
  double lut_courbure_clust[3 * 256];
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.515600;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.531300;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.546900;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.562500;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.578100;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.593800;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.609400;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.625000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.640600;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.656300;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.671900;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.687500;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.703100;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.718800;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.734400;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.750000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.765600;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.781300;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.796900;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.812500;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.828100;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.843800;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.859400;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.875000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.890600;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.906300;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.921900;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.937500;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.953100;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.968800;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.984400;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.015600;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.031300;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.046900;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.062500;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.078100;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.093800;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.109400;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.125000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.140600;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.156300;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.171900;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.187500;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.203100;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.218800;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.234400;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.250000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.265600;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.281300;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.296900;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.312500;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.328100;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.343800;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.359400;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.375000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.390600;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.406300;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.421900;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.437500;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.453100;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.468800;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.484400;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.500000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.515600;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.531300;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.546900;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.562500;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.578100;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.593800;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.609400;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.625000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.640600;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.656300;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.671900;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.687500;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.703100;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.718800;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.734400;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.750000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.765600;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.781300;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.796900;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.812500;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.828100;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.843800;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.859400;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.875000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.890600;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.906300;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.921900;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.937500;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.953100;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.968800;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.984400;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.015600;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.031300;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.984400;
  lut_courbure_clust[i++] = 0.046900;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.968800;
  lut_courbure_clust[i++] = 0.062500;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.953100;
  lut_courbure_clust[i++] = 0.078100;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.937500;
  lut_courbure_clust[i++] = 0.093800;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.921900;
  lut_courbure_clust[i++] = 0.109400;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.906300;
  lut_courbure_clust[i++] = 0.125000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.890600;
  lut_courbure_clust[i++] = 0.140600;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.875000;
  lut_courbure_clust[i++] = 0.156300;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.859400;
  lut_courbure_clust[i++] = 0.171900;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.843800;
  lut_courbure_clust[i++] = 0.187500;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.828100;
  lut_courbure_clust[i++] = 0.203100;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.812500;
  lut_courbure_clust[i++] = 0.218800;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.796900;
  lut_courbure_clust[i++] = 0.234400;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.781300;
  lut_courbure_clust[i++] = 0.250000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.765600;
  lut_courbure_clust[i++] = 0.265600;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.750000;
  lut_courbure_clust[i++] = 0.281300;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.734400;
  lut_courbure_clust[i++] = 0.296900;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.718800;
  lut_courbure_clust[i++] = 0.312500;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.703100;
  lut_courbure_clust[i++] = 0.328100;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.687500;
  lut_courbure_clust[i++] = 0.343800;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.671900;
  lut_courbure_clust[i++] = 0.359400;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.656300;
  lut_courbure_clust[i++] = 0.375000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.640600;
  lut_courbure_clust[i++] = 0.390600;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.625000;
  lut_courbure_clust[i++] = 0.406300;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.609400;
  lut_courbure_clust[i++] = 0.421900;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.593800;
  lut_courbure_clust[i++] = 0.437500;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.578100;
  lut_courbure_clust[i++] = 0.453100;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.562500;
  lut_courbure_clust[i++] = 0.468800;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.546900;
  lut_courbure_clust[i++] = 0.484400;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.531300;
  lut_courbure_clust[i++] = 0.500000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.515600;
  lut_courbure_clust[i++] = 0.515600;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.500000;
  lut_courbure_clust[i++] = 0.531300;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.484400;
  lut_courbure_clust[i++] = 0.546900;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.468800;
  lut_courbure_clust[i++] = 0.562500;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.453100;
  lut_courbure_clust[i++] = 0.578100;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.437500;
  lut_courbure_clust[i++] = 0.593800;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.421900;
  lut_courbure_clust[i++] = 0.609400;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.406300;
  lut_courbure_clust[i++] = 0.625000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.390600;
  lut_courbure_clust[i++] = 0.640600;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.375000;
  lut_courbure_clust[i++] = 0.656300;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.359400;
  lut_courbure_clust[i++] = 0.671900;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.343800;
  lut_courbure_clust[i++] = 0.687500;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.328100;
  lut_courbure_clust[i++] = 0.703100;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.312500;
  lut_courbure_clust[i++] = 0.718800;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.296900;
  lut_courbure_clust[i++] = 0.734400;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.281300;
  lut_courbure_clust[i++] = 0.750000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.265600;
  lut_courbure_clust[i++] = 0.765600;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.250000;
  lut_courbure_clust[i++] = 0.781300;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.234400;
  lut_courbure_clust[i++] = 0.796900;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.218800;
  lut_courbure_clust[i++] = 0.812500;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.203100;
  lut_courbure_clust[i++] = 0.828100;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.187500;
  lut_courbure_clust[i++] = 0.843800;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.171900;
  lut_courbure_clust[i++] = 0.859400;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.156300;
  lut_courbure_clust[i++] = 0.875000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.140600;
  lut_courbure_clust[i++] = 0.890600;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.125000;
  lut_courbure_clust[i++] = 0.906300;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.109400;
  lut_courbure_clust[i++] = 0.921900;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.093800;
  lut_courbure_clust[i++] = 0.937500;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.078100;
  lut_courbure_clust[i++] = 0.953100;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.062500;
  lut_courbure_clust[i++] = 0.968800;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.046900;
  lut_courbure_clust[i++] = 0.984400;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.031300;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.015600;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.984400;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.968800;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.953100;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.937500;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.921900;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.906300;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.890600;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.875000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.859400;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.843800;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.828100;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.812500;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.796900;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.781300;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.765600;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.750000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.734400;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.718800;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.703100;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.687500;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.671900;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.656300;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.640600;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.625000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.609400;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.593800;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.578100;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.562500;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.546900;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.531300;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.515600;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.500000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.484400;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.468800;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.453100;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.437500;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.421900;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.406300;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.390600;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.375000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.359400;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.343800;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.328100;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.312500;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.296900;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.281300;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.265600;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.250000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.234400;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.218800;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.203100;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.187500;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.171900;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.156300;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.140600;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.125000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.109400;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.093800;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.078100;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.062500;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.046900;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.031300;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.015600;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 1.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.984400;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.968800;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.953100;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.937500;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.921900;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.906300;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.890600;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.875000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.859400;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.843800;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.828100;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.812500;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.796900;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.781300;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.765600;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.750000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.734400;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.718800;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.703100;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.687500;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.671900;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.656300;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.640600;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.625000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.609400;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.593800;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.578100;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.562500;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.546900;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.531300;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.515600;
  lut_courbure_clust[i++] = 0.000000;
  lut_courbure_clust[i++] = 0.000000;


  // read mesh from file
  FEVV::MeshSurface m_original;
  FEVV::PMapsContainer pmaps_bag_original;
  FEVV::Filters::read_mesh(
      input_original_file_path, m_original, pmaps_bag_original);

  FEVV::MeshSurface m_degraded;
  FEVV::PMapsContainer pmaps_bag_degraded;
  FEVV::Filters::read_mesh(
      input_deformed_file_path, m_degraded, pmaps_bag_degraded);

  auto pm_degrad = get(boost::vertex_point, m_degraded);
  auto pm_original = get(boost::vertex_point, m_original);

  using FaceNormalMap =
      typename FEVV::PMap_traits< FEVV::face_normal_t,
                                  FEVV::MeshSurface >::pmap_type;
  FaceNormalMap fnm_degrad;
  if(has_map(pmaps_bag_degraded, FEVV::face_normal))
  {
    std::cout << "use existing face-normal map for degraded mesh" << std::endl;
    fnm_degrad =
        get_property_map(FEVV::face_normal, m_degraded, pmaps_bag_degraded);
  }
  else
  {
    std::cout << "create face-normal map for degraded mesh" << std::endl;
    fnm_degrad = make_property_map(FEVV::face_normal, m_degraded);
    // store property map in property maps bag
    put_property_map(
        FEVV::face_normal, m_degraded, pmaps_bag_degraded, fnm_degrad);
    FEVV::Filters::calculate_face_normals(m_degraded, pm_degrad, fnm_degrad);
  }

  FaceNormalMap fnm_original;
  if(has_map(pmaps_bag_original, FEVV::face_normal))
  {
    std::cout << "use existing face-normal map for original mesh" << std::endl;
    fnm_original =
        get_property_map(FEVV::face_normal, m_original, pmaps_bag_original);
  }
  else
  {
    std::cout << "create face-normal map for original mesh" << std::endl;
    fnm_original = make_property_map(FEVV::face_normal, m_original);
    // store property map in property maps bag
    put_property_map(
        FEVV::face_normal, m_original, pmaps_bag_original, fnm_original);
    FEVV::Filters::calculate_face_normals(
        m_original, pm_original, fnm_original);
  }


  typedef typename FEVV::Vertex_pmap< FEVV::MeshSurface, double >
      VertexMSDM2Map;
  VertexMSDM2Map msdm2_pmap;

  if(FEVV::has_map(pmaps_bag_degraded, std::string("v:msdm2")))
  {
    msdm2_pmap =
        boost::any_cast< VertexMSDM2Map >(pmaps_bag_degraded.at("v:msdm2"));
  }
  else
  {
    msdm2_pmap =
        FEVV::make_vertex_property_map< FEVV::MeshSurface, double >(m_degraded);
    pmaps_bag_degraded["v:msdm2"] = msdm2_pmap;
  }

  int nb_levels = 3;
  double msd_m2;

  FEVV::Filters::process_msdm2_multires(m_degraded,
                                        pm_degrad,
                                        fnm_degrad,
                                        pmaps_bag_degraded,
                                        m_original,
                                        pm_original,
                                        fnm_original,
                                        pmaps_bag_original,
                                        nb_levels,
                                        msdm2_pmap,
                                        msd_m2);

  std::cout << "Calculated MSDM2 : " << msd_m2 << std::endl;

  using VertexColorMap =
      typename FEVV::PMap_traits< FEVV::vertex_color_t,
                                  FEVV::MeshSurface >::pmap_type;
  VertexColorMap v_cm;
  if(has_map(pmaps_bag_degraded, FEVV::vertex_color))
  {
    std::cout << "use existing vertex-color map" << std::endl;
    v_cm = get_property_map(FEVV::vertex_color, m_degraded, pmaps_bag_degraded);
  }
  else
  {
    std::cout << "create vertex-color map" << std::endl;
    v_cm = make_property_map(FEVV::vertex_color, m_degraded);
    // store property map in property maps bag
    put_property_map(FEVV::vertex_color, m_degraded, pmaps_bag_degraded, v_cm);
  }

  double max_msd_m2, min_msd_m2;

  FEVV::Filters::compute_min_max_vertices(
      m_degraded, msdm2_pmap, min_msd_m2, max_msd_m2);

  FEVV::Filters::color_vertices_from_map(
      m_degraded, msdm2_pmap, v_cm, min_msd_m2, max_msd_m2, lut_courbure_clust);

  FEVV::Filters::write_mesh(output_file_path, m_degraded, pmaps_bag_degraded);

  // check output file
  std::cout << "Comparing output file '" << output_file_path
            << "' with reference file '" << reference_file_path << "'..."
            << std::endl;

  if(FEVV::FileUtils::has_extension(output_file_path, ".off") ||
     FEVV::FileUtils::has_extension(output_file_path, ".coff"))
  {
    // use OFF file comparator
    if(!are_meshes_equal(output_file_path, reference_file_path, false))
    {
      std::cout << "Files are different!" << std::endl;
      return EXIT_FAILURE;
    }
  }
  else
  {
    // use text file comparator
    if(!identical_text_based_files(output_file_path, reference_file_path))
    {
      std::cout << "Files are different!" << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::cout << "Files are identical." << std::endl;

  return 0;
}
