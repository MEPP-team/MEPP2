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

#include "FEVV/Filters/Generic/color_mesh.h" // Required to visualize CMDM as a color map
#include "FEVV/Filters/Generic/minmax_map.h" // Required to visualize CMDM as a color map

#include <chrono> //time computation

//---------------------------------------------------------
//                       Main
//---------------------------------------------------------

// Main: load a mesh, apply the filter, write the mesh
int
main(int narg, char **argv)
{
  // input and output files
  std::string input_original;
  std::string input_degraded;

  std::string output_file_path = "cmdm.output.obj";
  
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();//Chrono start

  if(narg == 3)
  {
    std::string reader = std::string(argv[1]);
    input_original = reader;
    std::string reader2 = std::string(argv[2]);
    input_degraded = reader2;
  }
  else
  {
    std::cout << "Use of " << argv[0] << " is :\n";
    std::cout << "\t" << argv[0]
              << " path/to/original_mesh path/to/degraded_mesh" << std::endl;
    std::cout << argv[0]
              << " will compute the CMDM value between the degraded mesh and "
                 "the original.\n"
              << std::endl;
    return 1;
  }

  // create LUT
  FEVV::Filters::ColorMeshLUT lut_courbure_clust = FEVV::Filters::make_LUT(false);

  // read mesh from file
  FEVV::MeshSurface m_original;
  FEVV::PMapsContainer pmaps_bag_original;
  FEVV::Filters::read_mesh(input_original, m_original, pmaps_bag_original);

  FEVV::MeshSurface m_degraded;
  FEVV::PMapsContainer pmaps_bag_degraded;
  FEVV::Filters::read_mesh(input_degraded, m_degraded, pmaps_bag_degraded);

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
    std::cout <<
           "CMDM needs 2 meshes with diffuse colors.\nUse MSDM2 Plugin for "
           "meshes having geometry characteristics only."
        << std::endl;
    return EXIT_FAILURE;
  }

  // Using a vertex property map to store the local distortions and visualize them later (color map).
  typedef typename FEVV::Vertex_pmap< FEVV::MeshSurface , double > vertex_cmdm_map;
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
  double  cmdm_1_2, cmdm_2_1, cmdm;

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

 FEVV::Filters::process_CMDM_multires(m_original,
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
 
  cmdm = (cmdm_1_2 + cmdm_2_1) / 2.;

  std::cout << "Calculated CMDM = " << cmdm << std::endl;
  
  // Visalize local distortions of the degraded mesh (color map of cmdm_1_2)
  using VertexColorMap =
      typename FEVV::PMap_traits< FEVV::vertex_color_t,
                                  FEVV::MeshSurface >::pmap_type;
  VertexColorMap v_cm =
      get_property_map(FEVV::vertex_color, m_degraded, pmaps_bag_degraded);

  double max_cmdm_1_2, min_cmdm_1_2;

  FEVV::Filters::compute_min_max_vertices(m_degraded, 
										  cmdm_pmap_deg, 
										  min_cmdm_1_2, 
										  max_cmdm_1_2);

  FEVV::Filters::color_vertices_from_map(m_degraded,
                                         cmdm_pmap_deg,
                                         v_cm,
                                         min_cmdm_1_2,
                                         max_cmdm_1_2,
                                         lut_courbure_clust);

  FEVV::Filters::write_mesh(output_file_path, m_degraded, pmaps_bag_degraded);

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now(); //Chrono end
  double time = std::chrono::duration_cast< std::chrono::microseconds > (end - begin).count() / 1000000.0;
  std::cout << time << " sec elapsed " << std::endl;
  
  return 0;
}
