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

#include "FEVV/Wrappings/Geometry_traits.h"


namespace FEVV {
namespace Filters {


/**
 * \brief  Normalize each vector of a property map, one by one.
 *         Generic implementation.
 *
 * \param  begin       Begin iterator of the key of the property map
 * \param  end         End iterator of the key of the property map
 * \param  prop_map    Property map containing vectors
 */
template< typename Iterator,
          typename PropertyMap,
          typename GeometryTraits >
void
normalize_vector_map(const Iterator       &begin,
                     const Iterator       &end,
                     PropertyMap          &prop_map,
                     const GeometryTraits &gt)
{
  for(Iterator it = begin; it != end; ++it)
  {
    auto vec = get(prop_map, *it);
    auto vec_normalized = gt.normalize(vec);
    put(prop_map, *it, vec_normalized);
  }
}


/**
 * \brief  Normalize each vector of the vertex property map.
 *
 * \param  g                 Mesh with the property map
 * \param  prop_map          Property map of vectors
 */
template< typename HalfedgeGraph,
          typename PropertyMap >
void
normalize_vector_map_vertices(const HalfedgeGraph &g,
                              PropertyMap         &prop_map)
{
  FEVV::Geometry_traits< HalfedgeGraph > gt(g);
  normalize_vector_map(vertices(g).first,
                       vertices(g).second,
                       prop_map,
                       gt);
}


/**
 * \brief  Normalize each vector of the face property map.
 *
 * \param  g                 Mesh with the property map
 * \param  prop_map          Property map of vectors
 */
template< typename HalfedgeGraph,
          typename PropertyMap >
void
normalize_vector_map_faces(const HalfedgeGraph &g,
                           PropertyMap         &prop_map)
{
  FEVV::Geometry_traits< HalfedgeGraph > gt(g);
  normalize_vector_map(faces(g).first,
                       faces(g).second,
                       prop_map,
                       gt);
}


/**
 * \brief  Normalize each vector of the edge property map.
 *
 * \param  g                 Mesh with the property map
 * \param  prop_map          Property map of vectors
 */
template< typename HalfedgeGraph,
          typename PropertyMap >
void
normalize_vector_map_edges(const HalfedgeGraph &g,
                           PropertyMap         &prop_map)
{
  FEVV::Geometry_traits< HalfedgeGraph > gt(g);
  normalize_vector_map(edges(g).first,
                       edges(g).second,
                       prop_map,
                       gt);
}


/**
 * \brief  Normalize each vector of the halfedge property map.
 *
 * \param  g                 Mesh with the property map
 * \param  prop_map          Property map of vectors
 */
template< typename HalfedgeGraph,
          typename PropertyMap >
void
normalize_vector_map_halfedges(const HalfedgeGraph &g,
                               PropertyMap         &prop_map)
{
  FEVV::Geometry_traits< HalfedgeGraph > gt(g);
  normalize_vector_map(halfedges(g).first,
                       halfedges(g).second,
                       prop_map,
                       gt);
}


} // namespace Filters
} // namespace FEVV
