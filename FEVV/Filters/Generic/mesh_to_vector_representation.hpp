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

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Wrappings/properties.h"

#include "FEVV/Types/Mesh_vector_representation.h"


namespace FEVV {
namespace Filters {


/**
 * \brief  Build the vector representation of the mesh.
 *
 * \param  g           mesh
 * \param  pmaps       property maps bag of the mesh
 * \param  mvr         the vector representation structure to fill (output)
 * \param  gt          the geometry traits to use
 *
 * \sa     the simplified variant that use the default geometry traits
 *         of the mesh.
 */
template< typename FaceListGraph,
          typename coordP_type,
          typename coordN_type,
          typename coordT_type,
          typename coordC_type,
          typename index_type,
          typename GeometryTraits >
void
mesh_to_vector_representation(const FaceListGraph  &g,
                              const FEVV::PMapsContainer &pmaps,
                              FEVV::Types::MVR < coordP_type,
                                                 coordN_type,
                                                 coordT_type,
                                                 coordC_type,
                                                 index_type > &mvr,
                              const GeometryTraits &/*gt*/)
{
  typedef boost::graph_traits< FaceListGraph >       GraphTraits;
  typedef typename GraphTraits::vertex_descriptor    vertex_descriptor;
  typedef typename GraphTraits::vertex_iterator      vertex_iterator;
  typedef boost::iterator_range< vertex_iterator >   vertex_range;
  //typedef typename GraphTraits::edge_descriptor      edge_descriptor;
  typedef typename GraphTraits::halfedge_descriptor  halfedge_descriptor;
  //typedef typename GraphTraits::face_descriptor      face_descriptor;
  typedef typename GraphTraits::face_iterator        face_iterator;
  typedef boost::iterator_range< face_iterator >     face_range;
  typedef typename GeometryTraits::Point             Point;
  typedef typename GeometryTraits::Vector            Vector;

  // which vertex attribute are available for export ?

  bool use_vertex_normal = false;
  bool use_vertex_color = false;
  bool use_vertex_texture_coord = false;
  bool use_halfedge_texture_coord = false;
  bool use_face_color = false;
  bool use_face_normal = false; // per-face vertices normals
  bool use_face_material = false;
  bool use_mesh_materials = false;

  if(has_map(pmaps, vertex_normal))
    use_vertex_normal = true;
  else if(has_map(pmaps, face_normal))
    use_face_normal = true;
  if(has_map(pmaps, vertex_color))
    use_vertex_color = true;
  if(has_map(pmaps, face_color))
    use_face_color = true;
  if(has_map(pmaps, vertex_texcoord))
    use_vertex_texture_coord = true;
  else if(has_map(pmaps, halfedge_texcoord))
    use_halfedge_texture_coord = true;
  if(has_map(pmaps, FEVV::face_material))
    use_face_material = true;
  if(has_map(pmaps, FEVV::mesh_materials))
    use_mesh_materials = true;

  long vertex_index = 0; // to be sure to start to 0
  std::map< vertex_descriptor, long >
      index_map; // for correct corresponding between vertex order in file and
                 // vertices around face indices

  // PUT VERTICES IN CONTAINER

  vertex_range v_range = vertices(g);
  auto point_pm = get(boost::vertex_point, g);

  // loop over vertices
  for(vertex_iterator it_v = v_range.begin(); it_v != v_range.end(); ++it_v)
  {
    Point p = get(point_pm, *it_v);
    std::vector< coordP_type > point;

    point.push_back(p[0]);
    point.push_back(p[1]);
    point.push_back(p[2]);

    mvr.points_coords.push_back(point);

    index_map[*it_v] = vertex_index;
    ++vertex_index;

    // deal with vertex attributes

    if(use_vertex_normal)
    {
      // store vertex normal in standard container
      auto vn_pm = get_property_map(FEVV::vertex_normal, g, pmaps);
      auto vn = get(vn_pm, *it_v);
      std::vector< coordN_type > normal;
      normal.push_back(vn[0]);
      normal.push_back(vn[1]);
      normal.push_back(vn[2]);
      mvr.normals_coords.push_back(normal);
    }

    if(use_vertex_color)
    {
      // store vertex color in standard container
      auto vc_pm = get_property_map(FEVV::vertex_color, g, pmaps);
      auto vc = get(vc_pm, *it_v);
      std::vector< coordC_type > color;
      color.push_back(static_cast< coordC_type >(vc[0]));
      color.push_back(static_cast< coordC_type >(vc[1]));
      color.push_back(static_cast< coordC_type >(vc[2]));
      mvr.vertex_color_coords.push_back(color);
    }

    if(use_vertex_texture_coord)
    {
      // store vertex texture-coordinates in standard container
      auto vt_pm = get_property_map(FEVV::vertex_texcoord, g, pmaps);
      auto vt = get(vt_pm, *it_v);
      // vt = (uvx, uvy, texture_index)
      // do NOT store texture index in output file
      std::vector< coordT_type > uvi;
      uvi.push_back(static_cast< coordT_type >(vt[0]));
      uvi.push_back(static_cast< coordT_type >(vt[1]));
      mvr.texture_coords.push_back(uvi);
    }
  }

  // PUT FACES IN CONTAINER

  // loop over faces
  face_range f_range = faces(g);
  for(face_iterator it_f = f_range.begin(); it_f != f_range.end(); ++it_f)
  {
    std::vector< index_type > face;
    std::vector< index_type > normal_indices;
    std::vector< index_type > texture_indices;

    // loop over vertices incident to face
    halfedge_descriptor h = halfedge(*it_f, g);
    vertex_descriptor vbegin = target(h, g);
    vertex_descriptor v = vbegin;
    unsigned int face_vertices_nbr = 0;
    do
    {
      face.push_back(
          index_map[v]); // for correct corresponding between vertex order in
                         // file and vertices around face indices
      face_vertices_nbr++;

      // face's vertices attributes (normal & texture)
      if(use_vertex_normal)
        normal_indices.push_back(index_map[v]);
      if(use_vertex_texture_coord)
        texture_indices.push_back(index_map[v]);
      else if(use_halfedge_texture_coord)
      {
        // retrieve halfedge texture-coordinates from property map
        auto ht_pm = get_property_map(FEVV::halfedge_texcoord, g, pmaps);
        auto ht = get(ht_pm, h);

        // store halfedge texture-coordinates in standard container
        // ht = (uvx, uvy, texture_index)
        // do NOT store texture index in output file
        std::vector< coordT_type > uvi;
        uvi.push_back(static_cast< coordT_type >(ht[0]));
        uvi.push_back(static_cast< coordT_type >(ht[1]));
        mvr.texture_coords.push_back(uvi);

        // store halfedge texture-coordinates index in standard container
        texture_indices.push_back(
            static_cast< index_type >(mvr.texture_coords.size() - 1));
      }

      // next vertex in the face
      h = next(h, g);
      v = target(h, g);
    } while(v != vbegin);

    // face attribute (normal)
    if(use_face_normal)
    {
      // retrieve face-normal from property map
      auto fn_pm = get_property_map(FEVV::face_normal, g, pmaps);
      auto fn = get(fn_pm, *it_f);

      // store face-normal in standard container
      std::vector< coordN_type > normal;
      normal.push_back(fn[0]);
      normal.push_back(fn[1]);
      normal.push_back(fn[2]);
      mvr.normals_coords.push_back(normal);

      // store normal index in standard container for each face vertex
      for(unsigned int i = 0; i < face_vertices_nbr; i++)
        normal_indices.push_back(
            static_cast< index_type >(mvr.normals_coords.size() - 1));
    }

    mvr.faces_indices.push_back(face);
    if(use_vertex_normal || use_face_normal)
      mvr.normal_face_indices.push_back(normal_indices);
    if(use_vertex_texture_coord || use_halfedge_texture_coord)
      mvr.texture_face_indices.push_back(texture_indices);

    // face color
    if(use_face_color)
    {
      // retrieve face-color from property map
      auto fc_pm = get_property_map(FEVV::face_color, g, pmaps);
      Vector fc = get(fc_pm, *it_f);

      // store face-color in standard container
      std::vector< coordC_type > color;
      color.push_back(static_cast< coordC_type >(fc[0]));
      color.push_back(static_cast< coordC_type >(fc[1]));
      color.push_back(static_cast< coordC_type >(fc[2]));
      mvr.face_color_coords.push_back(color);
    }

    // face material
    if(use_face_material)
    {
      // retrieve face-color from property map
      auto fm_pm = get_property_map(FEVV::face_material, g, pmaps);
      auto material_id = get(fm_pm, *it_f);

      // store face-material in standard container
      mvr.face_material.push_back(material_id);
    }
  }

  // PUT MATERIALS IN CONTAINER

  if(use_mesh_materials)
  {
    // loop over materials
    auto mm_pm = get_property_map(FEVV::mesh_materials, g, pmaps);
    auto it = mm_pm.storage_begin();
    auto it_end = mm_pm.storage_end();
    for(; it != it_end; ++it)
      mvr.materials.push_back(*it);
  }
}


/**
 * \brief  Build the vector representation of the mesh.
 *
 * \param  g           mesh
 * \param  pmaps       property maps bag of the mesh
 * \param  mvr         the vector representation structure to fill (output)
 *
 * \sa     the variant that use the geometry traits provided by the user.
 */
template< typename FaceListGraph,
          typename coordP_type,
          typename coordN_type,
          typename coordT_type,
          typename coordC_type,
          typename index_type,
          typename GeometryTraits = FEVV::Geometry_traits< FaceListGraph > >
void
mesh_to_vector_representation(const FaceListGraph  &g,
                              const FEVV::PMapsContainer &pmaps,
                              FEVV::Types::MVR < coordP_type,
                                                 coordN_type,
                                                 coordT_type,
                                                 coordC_type,
                                                 index_type > &mvr)
{
  GeometryTraits gt(g);
  mesh_to_vector_representation(g, pmaps, mvr, gt);
}


} // namespace Filters
} // namespace FEVV
