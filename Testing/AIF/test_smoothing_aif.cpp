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
#include "FEVV/Wrappings/Graph_traits_aif.h"
#include "FEVV/Wrappings/Geometry_traits_aif.h"

#include "FEVV/Wrappings/properties_aif.h"

using AIFKernel = FEVV::AIF_mesh_kernel_generator;
using MeshT = FEVV::DataStructures::AIF::AIFMesh;
using PtrMeshT = FEVV::DataStructures::AIF::AIFMesh::ptr_mesh;

//#define USE_GENERIC_READER_WRITER // when activated you need to link to vtk
//libs if VTK is used
#ifndef USE_GENERIC_READER_WRITER
#include "FEVV/DataStructures/AIF/AIFMeshReader.hpp"
#include "FEVV/DataStructures/AIF/AIFMeshWriter.hpp"
#else
//---------------- specific code above --------------------

//---------------- generic code below ---------------------

#include "FEVV/Filters/Generic/generic_reader.hpp"
#include "FEVV/Filters/Generic/generic_writer.hpp"
#endif

#include "FEVV/Filters/Generic/Manifold/calculate_vertices_one_ring_barycenter.hpp"
#include "FEVV/Filters/Generic/Manifold/calculate_vertices_one_ring_angles_based_centroid.hpp"
#include "FEVV/Filters/Generic/reposition_vertices.hpp"

//---------------- generic code above ---------------------

//---------------- other code below -----------------------

#include "FEVV/Tools/IO/FileUtilities.hpp"

#include <fstream>

void
test_calculate_scaling_aif(const std::string &input_file_path)
{
  using BarycenterMap =
      typename FEVV::Vertex_pmap_traits< MeshT, Point >::pmap_type;
  /**********************************************************************************************************/
#ifdef USE_GENERIC_READER_WRITER
  MeshT m;
  FEVV::PMapsContainer pmaps_bag;
  FEVV::Filters::read_mesh(input_file_path, m, pmaps_bag);
#else
  PtrMeshT pm;
  FEVV::DataStructures::AIF::AIFMeshReader in;
  if(!(pm = in.read(input_file_path)))
  {
    std::cout << "reading failed";
    return;
  }
  MeshT m(*pm); // copy contructor use
#endif
  /**********************************************************************************************************/
  try
  {
    /**********************************************************************************************************/
    auto pos_pm = get(boost::vertex_point, m);
    /**********************************************************************************************************/
    BarycenterMap barycenters_pm =
        FEVV::Vertex_pmap_traits< MeshT, Point >::create(m);
    FEVV::Filters::calculate_vertices_one_ring_barycenter(
        m, pos_pm, barycenters_pm, 0.2f);
    FEVV::Filters::reposition_vertices(m, pos_pm, barycenters_pm);

    FEVV::Filters::calculate_vertices_one_ring_angles_based_centroid(
        m, pos_pm, barycenters_pm, 0.2f);
    FEVV::Filters::reposition_vertices(m, pos_pm, barycenters_pm);
  }
  catch(const std::exception &e)
  {
    std::cout << "processing failed: " << e.what();
    return;
  }
#ifdef USE_GENERIC_READER_WRITER
  /**********************************************************************************************************/
  FEVV::Filters::write_mesh(
      "smoothed_aif_" + FEVV::FileUtils::get_file_name(input_file_path) +
          FEVV::FileUtils::get_file_extension(input_file_path),
      m,
      pmaps_bag);
#else
  /**********************************************************************************************************/
  FEVV::DataStructures::AIF::AIFMeshWriter out;
  try
  {
    out.write(m,
              "smoothed_aif_" +
                  FEVV::FileUtils::get_file_name(input_file_path) +
                  FEVV::FileUtils::get_file_extension(input_file_path));
  }
  catch(...)
  {
    std::cout << "writing failed";
    return;
  }
#endif
  std::cout << "Done." << std::endl;
}

//------------------------------------------------------------------------------

int
main(int narg, char **argv)
{
  if(narg < 2)
  {
    std::cout << "Usage: " << argv[0]
              << " filename; filename being an off file." << std::endl;
    exit(EXIT_FAILURE);
  }

  test_calculate_scaling_aif(argv[1]);
  return 0;
}
