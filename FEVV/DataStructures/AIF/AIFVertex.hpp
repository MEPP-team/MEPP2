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

#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include <limits>
#include <stdexcept>
#include <map>
#include <memory> // for unique_ptr
#include <array>
#include <boost/shared_ptr.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/property_map/property_map.hpp>


namespace FEVV {
namespace DataStructures {
namespace AIF {

class AIFEdge;
class AIFFace;

/**
 * \class	AIFVertex
 * \brief	This class represents a vertex used by AIFMesh objects.
 *			An AIFVertex natively saves relations with its incident
 *edges (AIFEdge). \see		AIFEdge \see		AIFFace
 */
class AIFVertex
{
public:
  friend class AIFTopologyHelpers;

public:
  typedef AIFVertex self;
  typedef boost::shared_ptr< self > ptr;
  typedef boost::shared_ptr< self > ptr_vertex;

  typedef AIFEdge edge_type;
  typedef AIFFace face_type;
  typedef boost::shared_ptr< edge_type > ptr_edge;
  typedef boost::shared_ptr< face_type > ptr_face;

  typedef std::vector< ptr_vertex > VertexContainerType;
  typedef std::vector< ptr_edge > EdgeContainerType;
  typedef std::vector< ptr_face > FaceContainerType;

  typedef double CoordinateType;
  typedef face_type::NormalCoordinateType NormalCoordinateType;

private:
  /*!
   *			AIFVertex unique index.
   */
  std::size_t m_Index;
  /*!
   *			Container on incident edges (AIFEdge pointer).
   */
  EdgeContainerType m_Incident_PtrEdges;

  static const unsigned int m_d = 3;

private:
  /*!
   *			Container on incident faces (AIFFace pointer).
   *			This container is filled only on demand and is modify if
   *necessary.
   */
  FaceContainerType m_Incident_PtrFaces; // is managed via AIFTopologyHelpers
  /*!
   *			Incident_PtrFaces cache state indicator.
   *			True if the cache is valid.
   */
  bool m_Incident_PtrFaces_Computed; // is managed via AIFTopologyHelpers
  /*!
   *			Container of one-ring vertices (AIFVertex pointer).
   *			This container is fill only on demand and is modified if
   *necessary. Handled with a cache system.
   */
  VertexContainerType m_One_Ring_Vertices; // is managed via AIFTopologyHelpers
  /*!
   *			One-ring vertices cache state indicator.
   *			True if the cache is valid.
   */
  bool m_Is_One_Ring_Vertices_Computed; // is managed via AIFTopologyHelpers

private:
  /*!
   *			AIFVertex default constructor.
   */
  AIFVertex()
      : m_Index(-1), m_Incident_PtrEdges(), m_Incident_PtrFaces(),
        m_Incident_PtrFaces_Computed(false), m_One_Ring_Vertices(),
        m_Is_One_Ring_Vertices_Computed(false)
  {
  }

  /*!
   * 			AIFVertex copy constructor
   *			Be careful, it create a new AIFVertex without the same native
   *incidence relations AND with the same index.
   */
  AIFVertex(const self &other)
      : m_Index(other.m_Index), // same index for the time being (this should be
                                // the case only in the Clone method)
        m_Incident_PtrEdges(),  // incidence relations are not copied for
                                // vertices
        m_Incident_PtrFaces(),  // incidence relations are not copied for
                                // vertices
        m_Incident_PtrFaces_Computed(false), m_One_Ring_Vertices(),
        m_Is_One_Ring_Vertices_Computed(false)
  {
  }

public:
  /*!
   * 			Static constructor method
   * \return	A pointer to a new instance of AIFVertex class.
   */
  static ptr_vertex New();
  /*!
   * 			Static constructor method
   * \param	other	The vertex to copy.
   * \return	A pointer to a new instance of AIFVertex class.
   */
  static ptr_vertex New(const self &other);
  /*!
   * 			Index getter
   * \return	The index of the involved vertex.
   */
  std::size_t GetIndex() const { return m_Index; }
  /*!
   * 			Index setter
   * \param	The index of the involved vertex.
   */
  void SetIndex(std::size_t idx) { m_Index = idx; }
  /*!
   * 			Native incidence relations getter
   * \return	An iterator range on the incident edges (AIFEdge pointer) of the
   * involved vertex.
   */
  boost::iterator_range< EdgeContainerType::const_iterator > GetIncidentEdges()
  {
    return boost::make_iterator_range(m_Incident_PtrEdges.cbegin(),
                                      m_Incident_PtrEdges.cend());
  }

  /*!
   * 			Degree getter
   * \return	The number of incident edges of the involved vertex.
   */
  unsigned int GetDegree()
  {
    return static_cast< unsigned int >(m_Incident_PtrEdges.size());
  }

  /*!
   * 			Print the textual description of the involved vertex on the
   * standard output.
   */
  void Print() const;

  /*!
   * 			Less operator (natural ordering)
   * \param	other	The second vertex to compare with the current one.
   * \return True if the current AIFVertex is strictly less than the other one.
   *         Else return false;
   */
  bool operator<(const self &other) const { return this < &other; }
  /*!
   * 			equality operator
   * \param	other	The second vertex to compare with the current one.
   * \return True if the current AIFVertex is equal to other one.
   *         Else return false;
   */
  bool operator==(const self &other) const { return this == &other; }
  /*!
   * 			difference operator
   * \param	other	The second vertex to compare with the current one.
   * \return True if the current AIFVertex is different from other one.
   *         Else return false;
   */
  bool operator!=(const self &other) const { return !(*this == other); }
  /*!
   * 			Less or equals operator
   * \param	other	The second vertex to compare with the current one.
   * \return True if the current AIFVertex is less or equal to other one.
   *         Else return false;
   */
  bool operator<=(const self &other) const
  {
    return (*this < other) || (*this == other);
  }
  /*!
   * 			Greater or equals operator
   * \param	other	The second vertex to compare with the current one.
   * \return True if the current AIFVertex is greater or equal to other one.
   *         Else return false;
   */
  bool operator>=(const self &other) const { return !(*this < other); }
  /*!
   * 			Greater operator
   * \param	other	The second vertex to compare with the current one.
   * \return True if the current AIFVertex is strictly greater to other one.
   *         Else return false;
   */
  bool operator>(const self &other) const { return (other < *this); }
};

} // namespace AIF
} // namespace DataStructures
} // namespace FEVV


#include "FEVV/DataStructures/AIF/AIFVertex.inl"
