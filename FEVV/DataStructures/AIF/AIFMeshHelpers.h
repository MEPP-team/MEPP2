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

#include "FEVV/DataStructures/AIF/AIFMesh.hpp"
#include "FEVV/Tools/Math/MatrixOperations.hpp"

namespace FEVV {
namespace DataStructures {
namespace AIF {

/**
 * \class	AIFMeshHelpers
 * \brief	This class is an helper class associated to the AIFMesh
 *structure. It hosts helper functions that do not rely only on topological
 *functions but also on properties (point coordinates for example). \see
 *AIFMesh
 */
class AIFMeshHelpers
{
public:
  typedef AIFPropertiesHelpers::CoordinateType CoordinateType;
  typedef AIFPropertiesHelpers::NormalCoordinateType NormalCoordinateType;
  typedef AIFPropertiesHelpers::Point Point;
  typedef AIFPropertiesHelpers::Vector Vector;

  typedef AIFTopologyHelpers::ptr_mesh ptr_mesh;
  typedef AIFTopologyHelpers::vertex_descriptor vertex_descriptor;
  typedef AIFTopologyHelpers::edge_descriptor edge_descriptor;
  typedef AIFTopologyHelpers::vertex_container_in_vertex
      vertex_container_in_vertex;
  typedef AIFTopologyHelpers::edge_container_in_vertex edge_container_in_vertex;
  typedef AIFTopologyHelpers::face_container_in_vertex face_container_in_vertex;
  typedef AIFTopologyHelpers::vertex_container_in_edge vertex_container_in_edge;
  typedef AIFTopologyHelpers::edge_container_in_edge edge_container_in_edge;
  typedef AIFTopologyHelpers::face_container_in_edge face_container_in_edge;
  typedef AIFTopologyHelpers::vertex_container_in_face vertex_container_in_face;
  typedef AIFTopologyHelpers::edge_container_in_face edge_container_in_face;
  typedef AIFTopologyHelpers::face_container_in_face face_container_in_face;
  typedef AIFTopologyHelpers::vertex_container vertex_container;
  typedef AIFTopologyHelpers::edge_container edge_container;
  typedef AIFTopologyHelpers::face_container face_container;

  /*!
   *
   *  from GHARIAL commit 24d25d4b145cb50087747183adb9c0894ab116fb (19 sept.
   * 2016)
   */
  template< typename GeometryTraits >
  static bool is_a_T_junction_vertex(vertex_descriptor v,
                                     ptr_mesh mesh,
                                     const GeometryTraits &gt)
  {
    if(AIFTopologyHelpers::is_surface_interior_vertex(v))
      return false;
    // all incident edges are interior edges, thus we are sure
    // there is not T-junction at this position

    ////////////////////////////////////////////////////////////////////////////
    vertex_container_in_vertex oneR =
        AIFTopologyHelpers::get_ordered_one_ring_of_adjacent_vertices(v);
    auto it = oneR.begin();
    auto ite = oneR.end();
    decltype(it) it2;
    if(it == ite)
      return false;
    ////////////////////////////////////////////////////////////////////////////
    // We must find 3 aligned points including current point (this) =>
    // T-junction detection which does not work for general crack detection for
    // 3D surfaces [in that case we can apply an orthogonal projection for this
    // point, but it is sensitive to numerical precision...]. Anyway, a crack
    // can also be considered as wanted surface holes...
    while(it != ite)
    {
      it2 = it;
      ++it2;
      while(it2 != ite)
      {
        Point &vPoint = AIFPropertiesHelpers::get_point(mesh, v);
        Point &itPoint = AIFPropertiesHelpers::get_point(mesh, *it);
        Point &it2Point = AIFPropertiesHelpers::get_point(mesh, *it2);

        std::vector< double > vVector(vPoint.cbegin(), vPoint.cend());
        std::vector< double > itVector(itPoint.cbegin(), itPoint.cend());
        std::vector< double > it2Vector(it2Point.cbegin(), it2Point.cend());

        if(Math::Vector::are_aligned(vVector, itVector, it2Vector) &&
           // must be on both sides of v vertex
           (gt.dot_product(Math::Vector::sub(vVector, itVector),
                           Math::Vector::sub(vVector, it2Vector)) < 0.))
        {
          if(AIFTopologyHelpers::is_surface_border_edge(
                 AIFTopologyHelpers::common_edge(v, *it)) &&
             AIFTopologyHelpers::is_surface_border_edge(
                 AIFTopologyHelpers::common_edge(v, *it2)))
          {
            if(AIFTopologyHelpers::are_adjacent(*it, *it2))
            {
              // simple T-junction case (2 triangles)
              return true;
            }

            vertex_container_in_vertex tmp =
                AIFTopologyHelpers::get_ordered_one_ring_of_adjacent_vertices(
                    *it2);
            auto tmpI = tmp.begin();
            auto tmpIe = tmp.end();

            while(tmpI != tmpIe)
            {
              if(**tmpI == *v)
              {
                ++tmpI;
                continue;
              }

              Point &tmpIPoint = AIFPropertiesHelpers::get_point(mesh, *tmpI);
              std::vector< double > tmpIVector(tmpIPoint.cbegin(),
                                               tmpIPoint.cend());

              if(Math::Vector::are_aligned(vVector, itVector, tmpIVector))
              {
                if(AIFTopologyHelpers::are_adjacent(*tmpI, *it) &&
                   AIFTopologyHelpers::is_surface_border_edge(
                       AIFTopologyHelpers::common_edge(*tmpI, *it)))
                {
                  // more complex T-junction case (3 triangles)
                  return true;
                }

                // more complex T-junction case (4 and more triangles)
                // are not taken into account yet
              }

              ++tmpI;
            }
          }
        }

        ++it2;
      }
      ++it;
    }
    ////////////////////////////////////////////////////////////////////////////
    return false; // not found any collinear triplet
  }

  /*!
   * 			A 2 manifold vertex one ring
   * \param	v	The involving vertex
   * \return	A boolean telling if the one-ring of v is 2-manifold, considering
   * both the topology and the vertex position.
   */
  template< typename GeometryTraits >
  static bool is_one_ring_2_manifold(vertex_descriptor v,
                                     ptr_mesh mesh,
                                     const GeometryTraits &gt)
  {
    if(is_a_T_junction_vertex(v, mesh, gt))
      return false;

    return AIFTopologyHelpers::is_one_ring_2_manifold(v);
  }
};


} // namespace AIF
} // namespace DataStructures
} // namespace FEVV

