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

#include "FEVV/Filters/CGAL/Surface_mesh/msdm2.h"
#include "FEVV/Filters/Generic/minmax_map.h"
#include "FEVV/Filters/Generic/color_mesh.h"
#include "FEVV/Filters/Generic/calculate_face_normals.hpp"

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

  std::string output_file_path = "msdm2.output.obj";

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
              << " path/to/original_mesh path/to/degraged_mesh" << std::endl;
    std::cout << argv[0]
              << " will compute the MSDM2 value between the degraded mesh and "
                 "the original.\n"
              << std::endl;
    return 1;
  }

  // create LUT
  FEVV::Filters::ColorMeshLUT lut_courbure_clust =
      FEVV::Filters::make_LUT(false);

  // read mesh from file
  FEVV::MeshSurface m_original;
  FEVV::PMapsContainer pmaps_bag_original;
  FEVV::Filters::read_mesh(input_original, m_original, pmaps_bag_original);

  FEVV::MeshSurface m_degraded;
  FEVV::PMapsContainer pmaps_bag_degraded;
  FEVV::Filters::read_mesh(input_degraded, m_degraded, pmaps_bag_degraded);

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
  double msdm2;

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
                                        msdm2);

  std::cout << "Calculated MSDM2 : " << msdm2 << std::endl;

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

  double max_msdm2, min_msdm2;

  FEVV::Filters::compute_min_max_vertices(
      m_degraded, msdm2_pmap, min_msdm2, max_msdm2);

  FEVV::Filters::color_vertices_from_map(
      m_degraded, msdm2_pmap, v_cm, min_msdm2, max_msdm2, lut_courbure_clust);

  FEVV::Filters::write_mesh(output_file_path, m_degraded, pmaps_bag_degraded);

  return 0;
}
