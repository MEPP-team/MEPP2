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

// Warning: do NOT include this file outside of AIFMesh.hpp

namespace FEVV {
namespace DataStructures {
namespace AIF {

/**
 * \class	AIFPropertiesHelpers
 * \brief	This class is an helper class associated to the AIFMesh
 *structure. AIFPropertiesHelpers implements all the basic geometric function
 *used to manipulate an AIFMesh point coordinates, and to compute normals or
 *edge length. Any such geometric operation on the AIF structure need to be
 *realized with this helper class. \see		AIFMesh
 */
class AIFPropertiesHelpers
{
public:
  typedef FEVV::DataStructures::AIF::AIFMesh AIFMesh;
  typedef AIFMesh::CoordinateType CoordinateType;
  typedef AIFMesh::NormalCoordinateType NormalCoordinateType;
  typedef AIFMesh::Point Point;
  typedef AIFMesh::Vector Vector;
  typedef AIFTopologyHelpers::smart_ptr_mesh smart_ptr_mesh;
  typedef AIFTopologyHelpers::ptr_mesh ptr_mesh;
  typedef AIFTopologyHelpers::ref_mesh ref_mesh;
  typedef AIFTopologyHelpers::vertex_descriptor vertex_descriptor;
  typedef AIFTopologyHelpers::edge_descriptor edge_descriptor;

  /*!
   *		  Set the vertex coordinates.
   * \param  m  mesh
   * \param  v  vertex descriptor
   * \param  x  x coordinate
   * \param  y  y coordinate
   * \param  z  z coordinate
   */
  static void set_point(const ptr_mesh m,
                        vertex_descriptor v,
                        CoordinateType x,
                        CoordinateType y,
                        CoordinateType z)
  {
    m->SetProperty< AIFVertex::ptr, Point >(
        "v:point", v->GetIndex(), Point(x, y, z));
  }

  /*!
   *		  Set the vertex coordinates.
   * \param  m  mesh
   * \param  v  vertex descriptor
   * \param  x  x coordinate
   * \param  y  y coordinate
   * \param  z  z coordinate
   */
  static void set_point(const smart_ptr_mesh m,
                        vertex_descriptor v,
                        CoordinateType x,
                        CoordinateType y,
                        CoordinateType z)
  {
    set_point(m.get(), v, x, y, z);
  }

  /*!
   *		  Set the vertex coordinates.
   * \param  m  mesh
   * \param  v  vertex descriptor
   * \param  x  x coordinate
   * \param  y  y coordinate
   * \param  z  z coordinate
   */
  static void set_point(/*const*/ ref_mesh m,
                        vertex_descriptor v,
                        CoordinateType x,
                        CoordinateType y,
                        CoordinateType z)
  {
    m.SetProperty< AIFVertex::ptr, Point >(
        "v:point", v->GetIndex(), Point(x, y, z));
  }

  /*!
   *		  Get the vertex coordinates.
   * \param  m  mesh
   * \param  v  vertex descriptor
   * \return  the point that contain the vertex coordinates
   */
  static Point &get_point(const ptr_mesh m, vertex_descriptor v)
  {
    return m->GetProperty< AIFVertex::ptr, Point >("v:point", v->GetIndex());
  }
  /*!
   *		  Get the vertex coordinates.
   * \param  m  mesh
   * \param  v  vertex descriptor
   * \return  the point that contain the vertex coordinates
   */
  static Point &get_point(const smart_ptr_mesh m, vertex_descriptor v)
  {
    return get_point(m.get(), v);
  }
  /*!
   *		  Get the vertex coordinates.
   * \param  m  mesh
   * \param  v  vertex descriptor
   * \return  the point that contain the vertex coordinates
   */
  static Point &get_point(/*const*/ ref_mesh m, vertex_descriptor v)
  {
    return m.GetProperty< AIFVertex::ptr, Point >("v:point", v->GetIndex());
  }

  /*!
   *        Compute the normal (non unitary) of the plane defined by the three
   *        points.
   */
  static Vector
  compute_normal(const Point &p1, const Point &p2, const Point &p3)
  {
    // calculate two vectors from the three first points
    Vector vector1(p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]);
    Vector vector2(p1[0] - p3[0], p1[1] - p3[1], p1[2] - p3[2]);

    // take the cross product of the two vectors to get
    NormalCoordinateType nx = static_cast< NormalCoordinateType >(
        vector1[1] * vector2[2] - vector1[2] * vector2[1]);
    NormalCoordinateType ny = static_cast< NormalCoordinateType >(
        vector1[2] * vector2[0] - vector1[0] * vector2[2]);
    NormalCoordinateType nz = static_cast< NormalCoordinateType >(
        vector1[0] * vector2[1] - vector1[1] * vector2[0]);

    Vector normal(nx, ny, nz);

    return normal;
  }

  /*!
   *        Compute the normal (non unitary) of the plane defined by the three
   *        vertices.
   */
  template< typename mesh >
  static Vector compute_normal(const mesh m,
                               const vertex_descriptor v1,
                               const vertex_descriptor v2,
                               const vertex_descriptor v3)
  {
    const Point &p1 = get_point(m, v1);
    const Point &p2 = get_point(m, v2);
    const Point &p3 = get_point(m, v3);
    return compute_normal(p1, p2, p3);
  }

  /*!
   *        Compute the unit normal of the plane defined by the three
   *        points.
   */
  static Vector
  compute_unit_normal(const Point &p1, const Point &p2, const Point &p3)
  {
    Vector normal = compute_normal(p1, p2, p3);
    double length = normal.length();

    if(length > 1e-8)
    {
      normal[0] /= (NormalCoordinateType)length;
      normal[1] /= (NormalCoordinateType)length;
      normal[2] /= (NormalCoordinateType)length;
    }

    return normal;
  }

  /*!
   *        Compute the unit normal of the plane defined by the three
   *        vertices.
   */
  template< typename mesh >
  static Vector compute_unit_normal(const mesh m,
                                    const vertex_descriptor v1,
                                    const vertex_descriptor v2,
                                    const vertex_descriptor v3)
  {
    const Point &p1 = get_point(m, v1);
    const Point &p2 = get_point(m, v2);
    const Point &p3 = get_point(m, v3);
    return compute_unit_normal(p1, p2, p3);
  }

  /*!
   *        Compute the length of an edge.
   */
  template< typename mesh >
  static CoordinateType length(const mesh m, const edge_descriptor e)
  {
    vertex_descriptor v1 = AIFTopologyHelpers::incident_vertices(e)[0];
    vertex_descriptor v2 = AIFTopologyHelpers::incident_vertices(e)[1];

    const Point &p1 = get_point(m, v1);
    const Point &p2 = get_point(m, v2);

    AIFVector< CoordinateType > v(p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]);
    return v.length();
  }
};


} // namespace AIF
} // namespace DataStructures
} // namespace FEVV
