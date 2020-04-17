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
#pragma once


#include "FEVV/Wrappings/properties.h"
#include "FEVV/Types/Mesh_vector_representation.h"
#include "FEVV/Filters/Generic/mesh_to_vector_representation.hpp"
#include "FEVV/Filters/Generic/mesh_from_vector_representation.hpp"


namespace FEVV {
namespace Filters {


/**
 * \brief  Copy a source mesh into a target mesh. Copy standard properties too.
 *
 * \param  g_s       the mesh to copy from
 * \param  pmaps_s   the source property maps bag
 * \param  g_t       the mesh to copy to
 * \param  pmaps_t   the target property maps bag
 * \param  gt_s      the geometry traits of source mesh
 * \param  gt_t      the geometry traits of target mesh
 *
 * \sa     the simplified variant that use the default geometry traits
 *         of the meshes.
 */
//TODO-elo: return v2v, e2e, h2h, f2f maps
template< typename FaceListGraphS,
          typename FaceListGraphT,
          typename GeometryTraitsS,
          typename GeometryTraitsT >
void
copy_graph(const FaceListGraphS  &g_s,
           const PMapsContainer  &pmaps_s,
           FaceListGraphT        &g_t,
           PMapsContainer        &pmaps_t,
           const GeometryTraitsS &gt_s,
           const GeometryTraitsT &gt_t)
{
  // build the vector representation of the source mesh

  // TODO-elo  templatize this 4 types ? Extract them from the mesh type ?
  typedef double coordP_type; // position coordinate type
  typedef double coordN_type; // normal coordinate type
  typedef float coordC_type;  // color coordinate type
  typedef float coordT_type;  // texture coordinate type
  typedef size_t index_type;

  FEVV::Types::MVR< coordP_type,
                    coordN_type,
                    coordT_type,
                    coordC_type,
                    index_type > mvr;

  mesh_to_vector_representation(g_s, pmaps_s, mvr, gt_s);

  // build the target mesh from the vector representation

  unsigned int duplicated_vertices_nbr = 0;
  bool use_corner_texture_coord = has_map(pmaps_s, halfedge_texcoord);

  mesh_from_vector_representation(g_t,
                                  pmaps_t,
                                  duplicated_vertices_nbr,
                                  mvr,
                                  use_corner_texture_coord,
                                  gt_t);

  // display some information

  std::cout << "copy_graph(): input mesh has "
            << size_of_vertices(g_s) << " vertices, "
            << size_of_edges(g_s) << " edges, "
            << size_of_faces(g_s) << " faces."
            << std::endl;
  std::cout << "copy_graph(): input mesh has property maps: [";
  for(auto &name: list_property_maps(pmaps_s))
    std::cout << " " << name;
  std::cout << " ]" << std::endl;

  std::cout << "copy_graph(): output mesh has "
            << size_of_vertices(g_t) << " vertices, "
            << size_of_edges(g_t) << " edges, "
            << size_of_faces(g_t) << " faces."
            << std::endl;
  if(duplicated_vertices_nbr > 0)
  {
    std::cout << "copy_graph(): "
              << duplicated_vertices_nbr << " vertices were duplicated."
              << std::endl;
  }
  std::cout << "copy_graph(): output mesh has property maps: [";
  for(auto &name: list_property_maps(pmaps_t))
    std::cout << " " << name;
  std::cout << " ]" << std::endl;
}


/**
 * \brief  Copy a source mesh into a target mesh. Copy standard properties too.
 *
 * \param  g_s       the mesh to copy from
 * \param  pmaps_s   the source property maps bag
 * \param  g_t       the mesh to copy to
 * \param  pmaps_t   the target property maps bag
 * \param  gt_s      the geometry traits of source mesh
 * \param  gt_t      the geometry traits of target mesh
 *
 * \sa     the variant that use the geometry traits provided by the user.
 */
template< typename FaceListGraphS,
          typename FaceListGraphT,
          typename GeometryTraitsS = FEVV::Geometry_traits< FaceListGraphS >,
          typename GeometryTraitsT = FEVV::Geometry_traits< FaceListGraphT > >
void
copy_graph(const FaceListGraphS  &g_s,
           const PMapsContainer  &pmap_s,
           FaceListGraphT        &g_t,
           PMapsContainer        &pmap_t)
{
  GeometryTraitsS gt_s(g_s);
  GeometryTraitsT gt_t(g_t);
  copy_graph(g_s, pmap_s, g_t, pmap_t, gt_s, gt_t);
}


} // namespace Filters
} // namespace FEVV
