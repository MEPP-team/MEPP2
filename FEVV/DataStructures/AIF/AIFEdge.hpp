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

#ifndef Q_MOC_RUN // See: https://bugreports.qt-project.org/browse/QTBUG-22829
#include <boost/property_map/property_map.hpp>
#endif

#include <iostream>
#include <vector>
#include <array>
#include <utility>
#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include <boost/range/iterator_range.hpp>

#include "FEVV/DataStructures/AIF/AIFVertex.hpp"


namespace FEVV {
namespace DataStructures {
namespace AIF {
class AIFEdge;
class AIFFace;
} // namespace AIF
} // namespace DataStructures
} // namespace FEVV

namespace FEVV {
namespace DataStructures {
namespace AIF {

/**
 * \class	AIFEdge
 * \brief	This class represents an edge used by AIFMesh objects.
 *			An AIFEdge natively saves relations with its incident
 *vertices (AIFVertex) and with its incident faces (AIFFace). \see
 *AIFVertex \see		AIFFace
 */
class AIFEdge
{
public:
  friend class AIFEdgeComparator;
  friend class AIFTopologyHelpers;

public:
  typedef AIFEdge self;
  typedef boost::shared_ptr< self > ptr;
  typedef boost::shared_ptr< self > ptr_edge;

  typedef AIFVertex vertex_type;
  typedef AIFFace face_type;
  typedef boost::shared_ptr< vertex_type > ptr_vertex;
  typedef boost::shared_ptr< face_type > ptr_face;

  typedef std::array< ptr_vertex, 2 > VertexContainerType;
  typedef std::vector< ptr_edge > EdgeContainerType;
  typedef std::vector< ptr_face > FaceContainerType;

private:
  /*!
   *			AIFEdge unique index.
   */
  std::size_t m_Index;
  /*!
   *			Container on incident vertices (AIFVertices pointer).
   */
  VertexContainerType m_incident_vertices;
  /*!
   *			Container on incident faces (AIFFace pointer).
   */
  FaceContainerType m_incident_PtrFaces;

private:
  /*!
   *			AIFEdge default constructor.
   */
  AIFEdge() : m_Index(-1), m_incident_vertices(), m_incident_PtrFaces() {}

  /*!
   * 			AIFEdge copy constructor
   *			Be careful, it create a new AIFEdge with the same native
   *incidence relations.
   */
  AIFEdge(const self &other)
      : m_Index(other.m_Index),
        m_incident_vertices(), // vertex incidence relations are not copied
        m_incident_PtrFaces()  // face incidence relations are not copied
  {
  }

public:
  /*!
   * 			Static constructor method
   * \return	A pointer to a new instance of AIFEdge class.
   */
  static ptr_edge New();
  /*!
   * 			Static constructor method
   * \param	other	The edge to copy.
   * \return	A pointer to a new instance of AIFEdge class.
   */
  static ptr_edge New(const self &other);
  /*!
   * 			Native incidence relations getter
   * \return	An iterator range on the incident vertices (AIFVertex pointer) of
   * the involved edge.
   */
  boost::iterator_range< VertexContainerType::const_iterator >
  GetIncidentVertices()
  {
    return boost::make_iterator_range(m_incident_vertices.cbegin(),
                                      m_incident_vertices.cend());
  }
  /*!
   * 			Native incidence relations getter
   * \return	An iterator range on the incident faces (AIFFace pointer) of the
   * involved edge.
   */
  boost::iterator_range< FaceContainerType::const_iterator >
  GetIncidentFaces() const
  {
    return boost::make_iterator_range(m_incident_PtrFaces.cbegin(),
                                      m_incident_PtrFaces.cend());
  }

  /*!
   * 			Native incidence relations getter
   * \return	An iterator range on the incident faces (AIFFace pointer) of the
   * involved edge.
   */
  boost::iterator_range< FaceContainerType::iterator > GetIncidentFaces()
  {
    return boost::make_iterator_range(m_incident_PtrFaces.begin(),
                                      m_incident_PtrFaces.end());
  }

  /*!
   * 			Degree getter
   * \return	The number of incident faces of the involved edge.
   */
  unsigned int GetDegree()
  {
    return static_cast< unsigned int >(m_incident_PtrFaces.size());
  }
  /*!
   * 			begin/first vertex getter
   * \return	The begin/first vertex of the involved edge.
   */
  ptr_vertex get_first_vertex() const { return m_incident_vertices[0]; }

  /*!
   * 			end/second vertex getter
   * \return	The end/second vertex of the involved edge.
   */
  ptr_vertex get_second_vertex() const { return m_incident_vertices[1]; }

  /*!
   * 			Index getter
   * \return	The index of the involved edge.
   */
  std::size_t GetIndex() const { return m_Index; }
  /*!
   * 			Index setter
   * \param	The index of the involved edge.
   */
  void SetIndex(std::size_t idx) { m_Index = idx; }
  /*!
   * 			First vertex setter
   * \param	v1	The new first vertex of the involved edge.
   */
  void set_first_vertex(ptr_vertex V1) { m_incident_vertices[0] = V1; }
  /*!
   * 			Second vertex setter
   * \param	v2	The new second vertex of the involved edge.
   */
  void set_second_vertex(ptr_vertex V2) { m_incident_vertices[1] = V2; }
  /*!
   *			Print the textual description of the involved edge on the
   *standard output.
   */
  void Print() const;

  /*!
   * 			Less operator (natural ordering)
   * \param	other	The second edge to compare with the current one.
   * \return True if the current AIFEdge is strictly less than the other one.
   *         Else return false;
   */
  bool operator<(const self &other) const;
  /*!
   * 			equality operator
   * \param	other	The second edge to compare with the current one.
   * \return True if the current AIFEdge is equal to other one.
   *         Else return false;
   */
  bool operator==(const self &other) const;
  /*!
   * 			difference operator
   * \param	other	The second edge to compare with the current one.
   * \return True if the current AIFEdge is different from other one.
   *         Else return false;
   */
  bool operator!=(const self &other) const { return !(*this == other); }
  /*!
   * 			Less or equals operator
   * \param	other	The second edge to compare with the current one.
   * \return True if the current AIFEdge is less or equal to other one.
   *         Else return false;
   */
  bool operator<=(const self &other) const
  {
    return (*this < other) || (*this == other);
  }
  /*!
   * 			Greater or equals operator
   * \param	other	The second edge to compare with the current one.
   * \return True if the current AIFEdge is greater or equal to other one.
   *         Else return false;
   */
  bool operator>=(const self &other) const { return !(*this < other); }
  /*!
   * 			Greater operator
   * \param	other	The second edge to compare with the current one.
   * \return True if the current AIFEdge is strictly greater to other one.
   *         Else return false;
   */
  bool operator>(const self &other) const { return (other < *this); }
};

} // namespace AIF
} // namespace DataStructures
} // namespace FEVV


#include "FEVV/DataStructures/AIF/AIFEdge.inl"
