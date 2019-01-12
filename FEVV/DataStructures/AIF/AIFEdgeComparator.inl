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

namespace FEVV {
namespace DataStructures {
namespace AIF {

inline
bool
AIFEdgeComparator::operator()(const AIFEdge &e1, const AIFEdge &e2)
{
  // compare the min vertex of each edge
  AIFEdge::ptr_vertex minVertexE1 = std::min< AIFEdge::ptr_vertex >(
      e1.m_incident_vertices[0], e1.m_incident_vertices[1]);
  AIFEdge::ptr_vertex minVertexE2 = std::min< AIFEdge::ptr_vertex >(
      e2.m_incident_vertices[0], e2.m_incident_vertices[1]);
  if(minVertexE1 != minVertexE2)
    return minVertexE1 < minVertexE2;

  // if min vertices are equal, then compare the max vertex of each edge
  AIFEdge::ptr_vertex maxVertexE1 = std::max< AIFEdge::ptr_vertex >(
      e1.m_incident_vertices[0], e1.m_incident_vertices[1]);
  AIFEdge::ptr_vertex maxVertexE2 = std::max< AIFEdge::ptr_vertex >(
      e2.m_incident_vertices[0], e2.m_incident_vertices[1]);
  return maxVertexE1 < maxVertexE2;
}

} // namespace AIF
} // namespace DataStructures
} // namespace FEVV
