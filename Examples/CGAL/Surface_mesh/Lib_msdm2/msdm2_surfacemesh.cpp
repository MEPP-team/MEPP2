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

#include "msdm2_surfacemesh.h"

#include "FEVV/Filters/CGAL/Surface_mesh/msdm2.h"
#include "FEVV/Filters/Generic/minmax_map.h"
#include "FEVV/Filters/Generic/color_mesh.h"
#include "FEVV/Filters/Generic/calculate_face_normals.hpp"


namespace msdm2 {

/**
 * SDM2 Filter library entry point
 */
void msdm2_surfacemesh(const MeshT &mesh_degraded,  /* input  */
                        const MeshT &mesh_original,  /* input  */
                        const int   nb_levels,       /* input  */
                        double      &msdm2_value,    /* output */
                        Msdm2MapT   &msdm2_pmap)     /* output */
{
  // sanity check, 1 <= nb_levels <= 3
  if(nb_levels < 1 || nb_levels > 3)
  {
    std::cout << "msdm2_precompiled: nb_levels must be in [1;3]. Aborting."
              << std::endl;
    return;
  }

  // retrieve Point maps
  auto pm_degraded = get(boost::vertex_point, mesh_degraded);
  auto pm_original = get(boost::vertex_point, mesh_original);

  using FaceNormalMap =
      typename FEVV::PMap_traits< FEVV::face_normal_t,
                                  MeshT                >::pmap_type;

  // create face normal map for degraded mesh
  FaceNormalMap fnm_degraded;
  std::cout << "create face-normal map for degraded mesh" << std::endl;//TODO-rm
  fnm_degraded = make_property_map(FEVV::face_normal, mesh_degraded);

  // calculate faces normals of degraded mesh
  FEVV::Filters::calculate_face_normals(mesh_degraded, pm_degraded,
                                        fnm_degraded);

  // create face normal map for original mesh
  FaceNormalMap fnm_original;
  std::cout << "create face-normal map for original mesh" << std::endl;//TODO-rm
  fnm_original = make_property_map(FEVV::face_normal, mesh_original);

  // calculate faces normals of original mesh
  FEVV::Filters::calculate_face_normals(mesh_original, pm_original,
                                        fnm_original);

  // create property maps bag used by msmd2 filter
  FEVV::PMapsContainer pmaps_bag_degraded;
  FEVV::PMapsContainer pmaps_bag_original;

  // call msdm2 filter
  FEVV::Filters::process_msdm2_multires(mesh_degraded,
                                        pm_degraded,
                                        fnm_degraded,
                                        pmaps_bag_degraded,
                                        mesh_original,
                                        pm_original,
                                        fnm_original,
                                        pmaps_bag_original,
                                        nb_levels,
                                        msdm2_pmap,
                                        msdm2_value);
}

} // namespace msdm2
