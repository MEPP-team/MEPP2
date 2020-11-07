// Copyright (c) 2012-2020 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published 
// by the Free Software Foundation; either version 3 of the License, 
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

#include "FEVV/Filters/CGAL/Surface_mesh/CMDM/cmdm.h"
#include "FEVV/Filters/Generic/calculate_face_normals.hpp"

#include "FEVV/Tools/IO/FileUtilities.hpp"
#include "Testing/Utils/utils_are_meshes_identical.hpp"
#include "Testing/Utils/utils_identical_text_based_files.hpp"

#include "FEVV/Filters/Generic/color_mesh.h" 
#include "FEVV/Filters/Generic/minmax_map.h"

int
main(int argc, const char **argv)
{
  if(argc != 5)
  {
    std::cout
        << "Load 2 meshes from input files, compute their CMDM, then compare the output file (the degraded mesh covered by the local distortion map) with the reference file."
        << std::endl;
    std::cout << "Usage: " << argv[0]
              << "  input_original_mesh_path  input_degraded_mesh_path Output_LD_degraded_mesh DegradedMesh_LD_OutputReferenceValues"
              << std::endl;
    std::cout << "Example: " << argv[0]
              << " Fish_ref.obj  Fish_simplified.obj  Fish_simplified_cmdm.out.obj  Fish_cmdm_reference.obj"
              << std::endl;
    return EXIT_FAILURE;
  }

  std::string input_original_file_path = argv[1];
  std::string input_deformed_file_path = argv[2];
  std::string output_file_path = argv[3];
  std::string reference_file_path = argv[4];

   // create LUT
  FEVV::Filters::ColorMeshLUT lut_courbure_clust = FEVV::Filters::make_LUT(false);

  // read meshes from file
  FEVV::MeshSurface m_original;
  FEVV::PMapsContainer pmaps_bag_original;
  FEVV::Filters::read_mesh(input_original_file_path, m_original, pmaps_bag_original);

  FEVV::MeshSurface m_degraded;
  FEVV::PMapsContainer pmaps_bag_degraded;
  FEVV::Filters::read_mesh(input_deformed_file_path, m_degraded, pmaps_bag_degraded);

  auto pm_degrad = get(boost::vertex_point, m_degraded);
  auto pm_original = get(boost::vertex_point, m_original);

  using FaceNormalMap = typename FEVV::PMap_traits< FEVV::face_normal_t, FEVV::MeshSurface >::pmap_type;
  FaceNormalMap fnm_degrad;
  if(has_map(pmaps_bag_degraded, FEVV::face_normal))
  {
    std::cout << "use existing face-normal map for degraded mesh" << std::endl;
    fnm_degrad = get_property_map(FEVV::face_normal, m_degraded, pmaps_bag_degraded);
  }
  else
  {
    std::cout << "create face-normal map for degraded mesh" << std::endl;
    fnm_degrad = make_property_map(FEVV::face_normal, m_degraded);
    // store property map in property maps bag
    put_property_map(FEVV::face_normal, m_degraded, pmaps_bag_degraded, fnm_degrad);
    FEVV::Filters::calculate_face_normals(m_degraded, pm_degrad, fnm_degrad);
  }

  FaceNormalMap fnm_original;
  if(has_map(pmaps_bag_original, FEVV::face_normal))
  {
    std::cout << "use existing face-normal map for original mesh" << std::endl;
    fnm_original = get_property_map(FEVV::face_normal, m_original, pmaps_bag_original);
  }
  else
  {
    std::cout << "create face-normal map for original mesh" << std::endl;
    fnm_original = make_property_map(FEVV::face_normal, m_original);
    // store property map in property maps bag
    put_property_map(FEVV::face_normal, m_original, pmaps_bag_original, fnm_original);
    FEVV::Filters::calculate_face_normals(m_original, pm_original, fnm_original);
  }

  // Verify if the original and the degraded meshes have diffuse color maps
  if(has_map(pmaps_bag_original, FEVV::vertex_color) &&
     has_map(pmaps_bag_degraded, FEVV::vertex_color))
  {
    std::cout << "use existing color maps for original and degraded meshes"
              << std::endl;
  }
  else
  {
    std::cout
        << "CMDM needs 2 meshes with diffuse colors.\nUse MSDM2 Plugin for "
           "meshes having geometry characteristics only."
        << std::endl;
    return EXIT_FAILURE;
  }

  // Using a vertex property map to store the local distortions and visualize them later (color map).
  typedef typename FEVV::Vertex_pmap< FEVV::MeshSurface, double > vertex_cmdm_map;
  vertex_cmdm_map cmdm_pmap_deg, cmdm_pmap_orig;

  if(FEVV::has_map(pmaps_bag_degraded, std::string("v:cmdm")))
  {
    cmdm_pmap_deg =
        boost::any_cast< vertex_cmdm_map >(pmaps_bag_degraded.at("v:cmdm"));
  }
  else
  {
    cmdm_pmap_deg =
        FEVV::make_vertex_property_map< FEVV::MeshSurface, double >(m_degraded);
    pmaps_bag_degraded["v:cmdm"] = cmdm_pmap_deg;
  }

  if(FEVV::has_map(pmaps_bag_original, std::string("v:cmdm")))
  {
    cmdm_pmap_orig =
        boost::any_cast< vertex_cmdm_map >(pmaps_bag_original.at("v:cmdm"));
  }
  else
  {
    cmdm_pmap_orig =
        FEVV::make_vertex_property_map< FEVV::MeshSurface, double >(m_original);
    pmaps_bag_original["v:cmdm"] = cmdm_pmap_orig;
  }

  int nb_levels = 3;
  double cmdm_1_2, cmdm;//cmdm_2_1,

  FEVV::Filters::process_CMDM_multires(m_degraded,
                                       pm_degrad,
                                       fnm_degrad,
                                       pmaps_bag_degraded,
                                       m_original,
                                       pm_original,
                                       fnm_original,
                                       pmaps_bag_original,
                                       nb_levels,
                                       cmdm_pmap_deg,
                                       cmdm_1_2);

  /*FEVV::Filters::process_CMDM_multires(m_original,
                                       pm_original,
                                       fnm_original,
                                       pmaps_bag_original,
                                       m_degraded,
                                       pm_degrad,
                                       fnm_degrad,
                                       pmaps_bag_degraded,
                                       nb_levels,
                                       cmdm_pmap_orig,
                                       cmdm_2_1);

  cmdm = (cmdm_1_2 + cmdm_2_1) / 2.;*/
  cmdm = cmdm_1_2;
  std::cout << "Calculated CMDM = " << cmdm << std::endl;

  // Visalize local distortions of the degraded mesh (color map of cmdm_1_2)
  using VertexColorMap =
      typename FEVV::PMap_traits< FEVV::vertex_color_t,
                                  FEVV::MeshSurface >::pmap_type;
  VertexColorMap v_cm =
      get_property_map(FEVV::vertex_color, m_degraded, pmaps_bag_degraded);

  double max_cmdm_1_2, min_cmdm_1_2;

  FEVV::Filters::compute_min_max_vertices(
      m_degraded, cmdm_pmap_deg, min_cmdm_1_2, max_cmdm_1_2);

  FEVV::Filters::color_vertices_from_map(m_degraded,
                                         cmdm_pmap_deg,
                                         v_cm,
                                         min_cmdm_1_2,
                                         max_cmdm_1_2,
                                         lut_courbure_clust);
  
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
