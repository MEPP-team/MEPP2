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

#include <map>


namespace FEVV {
namespace Filters {


/**
 * Helper type for vertex to vertex map.
 */
template< typename BoostGraphS,
          typename BoostGraphT >
using VertexToVertexMap = std::map<
   typename boost::graph_traits< BoostGraphS >::vertex_descriptor,
   typename boost::graph_traits< BoostGraphT >::vertex_descriptor >;


/**
 * Helper type for face to face map.
 */
template< typename FaceGraphS,
          typename FaceGraphT >
using FaceToFaceMap = std::map<
   typename boost::graph_traits< FaceGraphS >::face_descriptor,
   typename boost::graph_traits< FaceGraphT >::face_descriptor >;


/**
 * Implement the named parameters idiom for copy_graph()
 *
 * \param  v2v  vertex to vertex map to be populated
 * \param  f2f  face to face map to be populated
 */
template< typename FaceGraphS,
          typename FaceGraphT >
struct CopyGraphParameters
{
  typedef  VertexToVertexMap< FaceGraphS, FaceGraphT >  V2VMap;
  typedef  FaceToFaceMap< FaceGraphS, FaceGraphT >      F2FMap;

  CopyGraphParameters& v2v(V2VMap *_v2v)
  {
    m_v2v = _v2v;
    return *this;
  }

  CopyGraphParameters& f2f(F2FMap *_f2f)
  {
    m_f2f = _f2f;
    return *this;
  }

  V2VMap *m_v2v = nullptr;
  F2FMap *m_f2f = nullptr;
};


/**
 * \brief  Copy a source mesh into a target mesh. Copy standard properties too.
 *
 * \param  g_s       the mesh to copy from
 * \param  pmaps_s   the source property maps bag
 * \param  g_t       the mesh to copy to
 * \param  pmaps_t   the target property maps bag
 * \param  params    the named parameters
 * \param  gt_s      the geometry traits of source mesh
 * \param  gt_t      the geometry traits of target mesh
 *
 * \sa     the simplified variant that use the default geometry traits
 *         of the meshes.
 */
//TODO-elo: return v2v, e2e, h2h, f2f maps
template< typename FaceListGraphS,
          typename FaceListGraphT,
          typename Parameters,
          typename GeometryTraitsS,
          typename GeometryTraitsT >
void
copy_graph(const FaceListGraphS  &g_s,
           const PMapsContainer  &pmaps_s,
           FaceListGraphT        &g_t,
           PMapsContainer        &pmaps_t,
           const Parameters      &params,
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
  VertDescVect< FaceListGraphT > vd_target;
  FaceDescVect< FaceListGraphT > fd_target;
  bool use_corner_texture_coord = has_map(pmaps_s, halfedge_texcoord);

  MeshFromVectorReprParameters< FaceListGraphT > mfv_params;

  mesh_from_vector_representation(g_t,
                                  pmaps_t,
                                  duplicated_vertices_nbr,
                                  mvr,
                                  mfv_params.vd_target(&vd_target)
                                            .fd_target(&fd_target)
                                            .use_corner_texcoord(
                                                use_corner_texture_coord),
                                  gt_t);
 
  // populate v2v and f2f maps if provided
  
  // populate v2v map
  if(params.m_v2v)
  {
    // reset v2v map
    params.m_v2v->clear();

    // ensure there are as many source vertices than target vertices
    if(size_of_vertices(g_s) == vd_target.size())
    {
      // loop over source mesh vertices
      auto v_iter_pair = vertices(g_s);
      auto vi          = v_iter_pair.first;
      auto vi_end      = v_iter_pair.second;
      size_t i = 0;
      for(; vi != vi_end; ++vi)
      {
        // source vertices and vd_target are in the same order
        (*params.m_v2v)[*vi] = vd_target[i];
        i++;
      }
    }
    else
    {
      std::cout << "copy_graph(): WARNING unable to build v2v map because "
                   " source and target have not the same number of vertices."
                << std::endl;
    }
  }

  // populate f2f map
  if(params.m_f2f)
  {
    // reset f2f map
    params.m_f2f->clear();

    // ensure there are as many source faces than target faces
    if(size_of_faces(g_s) == fd_target.size())
    {
      // loop over source mesh faces
      auto f_iter_pair = faces(g_s);
      auto fi          = f_iter_pair.first;
      auto fi_end      = f_iter_pair.second;
      size_t i = 0;
      for(; fi != fi_end; ++fi)
      {
        // source faces and fd_target are in the same order
        (*params.m_f2f)[*fi] = fd_target[i];
        i++;
      }
    }
    else
    {
      std::cout << "copy_graph(): WARNING unable to build f2f map because "
                   " source and target have not the same number of faces."
                << std::endl;
    }
  }

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
 * \param  params    the named parameters
 *
 * \sa     the variant that use the geometry traits provided by the user.
 */
template< typename FaceListGraphS,
          typename FaceListGraphT,
          typename Parameters =
              CopyGraphParameters< FaceListGraphS, FaceListGraphT >,
          typename GeometryTraitsS = FEVV::Geometry_traits< FaceListGraphS >,
          typename GeometryTraitsT = FEVV::Geometry_traits< FaceListGraphT > >
void
copy_graph(const FaceListGraphS  &g_s,
           const PMapsContainer  &pmap_s,
           FaceListGraphT        &g_t,
           PMapsContainer        &pmap_t,
           const Parameters      &params =
               CopyGraphParameters< FaceListGraphS, FaceListGraphT >())
{
  GeometryTraitsS gt_s(g_s);
  GeometryTraitsT gt_t(g_t);
  copy_graph(g_s, pmap_s, g_t, pmap_t, params, gt_s, gt_t);
}


/**
 * \brief  Copy a source mesh into a target mesh. Do NOT Copy properties.
 *
 * \param  g_s       the mesh to copy from
 * \param  g_t       the mesh to copy to
 * \param  params    the named parameters
 * \param  gt_s      the geometry traits of source mesh
 * \param  gt_t      the geometry traits of target mesh
 *
 * \sa     the simplified variant that use the default geometry traits
 *         of the meshes.
 */
template< typename FaceListGraphS,
          typename FaceListGraphT,
          typename Parameters,
          typename GeometryTraitsS,
          typename GeometryTraitsT >
void
copy_graph(const FaceListGraphS  &g_s,
           FaceListGraphT        &g_t,
           const Parameters      &params,
           const GeometryTraitsS &gt_s,
           const GeometryTraitsT &gt_t)
{
  PMapsContainer  pmaps_s; // empty bag
  PMapsContainer  pmaps_t; // empty bag

  copy_graph(g_s, pmaps_s, g_t, pmaps_t, params, gt_s, gt_t);
}


/**
 * \brief  Copy a source mesh into a target mesh. Do NOT copy properties.
 *
 * \param  g_s       the mesh to copy from
 * \param  g_t       the mesh to copy to
 * \param  params    the named parameters
 *
 * \sa     the variant that use the geometry traits provided by the user.
 */
template< typename FaceListGraphS,
          typename FaceListGraphT,
          typename Parameters =
              CopyGraphParameters< FaceListGraphS, FaceListGraphT >,
          typename GeometryTraitsS = FEVV::Geometry_traits< FaceListGraphS >,
          typename GeometryTraitsT = FEVV::Geometry_traits< FaceListGraphT > >
void
copy_graph(const FaceListGraphS  &g_s,
           FaceListGraphT        &g_t,
           const Parameters      &params =
               CopyGraphParameters< FaceListGraphS, FaceListGraphT >())
{
  GeometryTraitsS gt_s(g_s);
  GeometryTraitsT gt_t(g_t);
  copy_graph(g_s, g_t, params, gt_s, gt_t);
}


} // namespace Filters
} // namespace FEVV
