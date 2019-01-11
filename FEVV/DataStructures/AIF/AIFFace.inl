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
#include "FEVV/DataStructures/AIF/AIFEdge.hpp"

namespace FEVV {
namespace DataStructures {
namespace AIF {

inline
AIFFace::ptr_face
AIFFace::New()
{
  ptr_face ptr(new self());
  return ptr;
}

inline
AIFFace::ptr_face
AIFFace::New(const self &other)
{
  ptr_face ptr(new self(other));
  return ptr;
}

inline
void
AIFFace::Print() const
{
  std::cout << "face " << this << " at index " << m_Index << std::endl;
  EdgeContainerType::const_iterator itBegin = m_Incident_PtrEdges.begin(),
                                    itEnd = m_Incident_PtrEdges.end();
  for(; itBegin != itEnd; ++itBegin)
    (*itBegin)->Print();
  std::cout << std::endl;
}

inline
boost::iterator_range< AIFFace::VertexContainerType::const_iterator >
AIFFace::GetIncidentVertices()
{
  if(m_Is_Incident_PtrVertices_Computed)
  {
    // adjacent vertices info is up-to-date, so we use it!
    // DBG std::cout << "Use face-incident-vertices cache for face " << this <<
    // std::endl;
    return m_Incident_PtrVertices;
  }

  // DBG std::cout << "No face-incident-vertices cache, computing it for face "
  // << this << std::endl;
  m_Incident_PtrVertices.clear();
  m_Is_Incident_PtrVertices_Computed = true;
  // we are computing a new one, so at the end of the computation this will be
  // true

  size_t nbEdge = m_Incident_PtrEdges.size();
  if(nbEdge == 1)
  { // Special case of a one edged face (only on face construction as this is an
    // invalid configuration)
    ptr_edge edge = *(m_Incident_PtrEdges.begin());
    m_Incident_PtrVertices.insert(m_Incident_PtrVertices.end(),
                                  edge->get_first_vertex());
    m_Incident_PtrVertices.insert(m_Incident_PtrVertices.end(),
                                  edge->get_second_vertex());
  }
  else if(nbEdge > 1)
  { // All other cases
    ptr_vertex currentVertex;
    EdgeContainerType::const_iterator itE = m_Incident_PtrEdges.cbegin();
    ptr_edge currentEdge = *itE;
    ++itE;
    ptr_edge nextEdge = *itE;

    if(currentEdge->get_first_vertex() == nextEdge->get_first_vertex() ||
       currentEdge->get_first_vertex() == nextEdge->get_second_vertex())
      currentVertex = currentEdge->get_first_vertex();
    else
      currentVertex = currentEdge->get_second_vertex();

    m_Incident_PtrVertices.insert(m_Incident_PtrVertices.end(), currentVertex);
    for(; itE != m_Incident_PtrEdges.cend(); ++itE)
    {
      if((*itE)->get_first_vertex() == currentVertex)
      {
        m_Incident_PtrVertices.insert(m_Incident_PtrVertices.end(),
                                      (*itE)->get_second_vertex());
        currentVertex = (*itE)->get_second_vertex();
      }
      else
      {
        m_Incident_PtrVertices.insert(m_Incident_PtrVertices.end(),
                                      (*itE)->get_first_vertex());
        currentVertex = (*itE)->get_first_vertex();
      }
    }
  } // if not the first two cases -> vertices container has been invalidate
    // after topological changes
  return boost::make_iterator_range(m_Incident_PtrVertices.cbegin(),
                                    m_Incident_PtrVertices.cend());
}

} // namespace AIF
} // namespace DataStructures
} // namespace FEVV
