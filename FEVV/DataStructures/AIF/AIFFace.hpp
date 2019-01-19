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

#include <iostream>
#include <vector>
//#include <set>
#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include <boost/range/iterator_range.hpp>


namespace FEVV {
namespace DataStructures {
namespace AIF {

class AIFVertex;
class AIFEdge;
class AIFFace;

} // namespace AIF
} // namespace DataStructures
} // namespace FEVV


namespace FEVV {
namespace DataStructures {
namespace AIF {

/**
 * \class	AIFFace
 * \brief	This class represents a face used by AIFMesh objects.
 *			An AIFFace natively saves relations with its incident edges
 *(AIFEdge). In addition, relations with its incident vertices (AIFVertex) are
 *saved on demand. \see		AIFVertex \see		AIFEdge
 */
class AIFFace
{
public:
  friend class AIFTopologyHelpers;

public:
  typedef AIFFace self;
  typedef boost::shared_ptr< self > ptr;
  typedef boost::shared_ptr< self > ptr_face;

  typedef AIFVertex vertex_type;
  typedef AIFEdge edge_type;
  typedef boost::shared_ptr< vertex_type > ptr_vertex;
  typedef boost::shared_ptr< edge_type > ptr_edge;

  typedef std::vector< ptr_vertex > VertexContainerType;
  typedef std::vector< ptr_edge > EdgeContainerType;
  typedef std::vector< ptr_face > FaceContainerType;

  typedef double NormalCoordinateType;

private:
  /*!
   *			AIFFace unique index.
   */
  std::size_t m_Index;
  /*!
   *			Container on incident edges (AIFEdge pointer).
   */
  EdgeContainerType m_Incident_PtrEdges;

private:
  /*!
   *			Container on incident vertices (AIFVertex pointer).
   *			This container is filled only on demand and is modified if
   *necessary.
   */
  VertexContainerType m_Incident_PtrVertices;

  /*!
   *			Face-incident-vertices cache validity flag.
   */
  bool m_Is_Incident_PtrVertices_Computed;

private:
  /*!
   *			AIFFace default constructor.
   */
  AIFFace()
      : m_Index(-1), m_Incident_PtrEdges(), m_Incident_PtrVertices(),
        m_Is_Incident_PtrVertices_Computed(false)
  {
  }
  /*!
   * 			AIFFace copy constructor
   *			Be careful, it create a new AIFFace with the same native
   *incidence relations AND the same index.
   */
  AIFFace(const self &other)
      : m_Index(other.m_Index),
        m_Incident_PtrEdges(),    // incidence relations are not copied
        m_Incident_PtrVertices(), // incidence relations are not copied
        m_Is_Incident_PtrVertices_Computed(false)
  {
  }

public:
  /*!
   * 			Static constructor method
   * \return	A pointer to a new instance of AIFFace class.
   */
  static ptr_face New();

  /*!
   * 			Static constructor method
   * \param	other	The face to copy.
   * \return	A pointer to a new instance of AIFFace class.
   */
  static ptr_face New(const self &other);

  /*!
   * 			Index getter
   * \return	The index of the involved face.
   */

  std::size_t GetIndex() const { return m_Index; }
  /*!
   * 			Index setter
   * \param	The index of the involved face.
   */
  void SetIndex(std::size_t idx) { m_Index = idx; }

  /*!
   * 			Native incidence relations getter
   * \return	An iterator range on the incident vertices (AIFVertex pointer) of
   * the involved face.
   */
  boost::iterator_range< VertexContainerType::const_iterator >
  GetIncidentVertices(); // cannot be const since it updates vertex incidency
                         // caching when necessary

  /*!
   * 			Native incidence relations getter
   * \return	An iterator range on the incident edges (AIFEdge pointer) of the
   * involved face.
   */
  boost::iterator_range< EdgeContainerType::const_iterator >
  GetIncidentEdges() const
  {
    return boost::make_iterator_range(m_Incident_PtrEdges.cbegin(),
                                      m_Incident_PtrEdges.cend());
  }

  /*!
   * 			Native incidence relations getter
   * \return	An iterator range on the incident edges (AIFEdge pointer) of the
   * involved face.
   */
  boost::iterator_range< EdgeContainerType::iterator > GetIncidentEdges()
  {
    return boost::make_iterator_range(m_Incident_PtrEdges.begin(),
                                      m_Incident_PtrEdges.end());
  }

  /*!
   * 			Degree getter
   * \return	The number of incident edges of the involved face.
   */
  unsigned int GetDegree()
  {
    return static_cast< unsigned int >(m_Incident_PtrEdges.size());
  }
  /*!
   * 			Print the textual description of the involved face on the standard
   * output.
   */
  void Print() const;

  /*!
   * 			Clear vertex incidency caching.
   */
  void clear_vertex_incidency()
  {
    // DBG std::cout << "clear face-incident-vertices cache for face " << this
    // << std::endl;
    m_Is_Incident_PtrVertices_Computed = false;
    m_Incident_PtrVertices.clear();
  }

  /*!
   * 			Less operator (natural ordering)
   * \param	other	The second face to compare with the current one.
   * \return True if the current AIFFace is strictly less than the other one.
   *         Else return false;
   */
  bool operator<(const self &other) const { return this < &other; }
  /*!
   * 			equality operator
   * \param	other	The second face to compare with the current one.
   * \return True if the current AIFFace is equal to other one.
   *         Else return false;
   */
  bool operator==(const self &other) const { return this == &other; }
  /*!
   * 			difference operator
   * \param	other	The second face to compare with the current one.
   * \return True if the current AIFFace is different from other one.
   *         Else return false;
   */
  bool operator!=(const self &other) const { return !(*this == other); }
  /*!
   * 			Less or equals operator
   * \param	other	The second face to compare with the current one.
   * \return True if the current AIFFace is less or equal to other one.
   *         Else return false;
   */
  bool operator<=(const self &other) const
  {
    return (*this < other) || (*this == other);
  }
  /*!
   * 			Greater or equals operator
   * \param	other	The second face to compare with the current one.
   * \return True if the current AIFFace is greater or equal to other one.
   *         Else return false;
   */
  bool operator>=(const self &other) const { return !(*this < other); }
  /*!
   * 			Greater operator
   * \param	other	The second face to compare with the current one.
   * \return True if the current AIFFace is strictly greater to other one.
   *         Else return false;
   */
  bool operator>(const self &other) const { return (other < *this); }
};

} // namespace AIF
} // namespace DataStructures
} // namespace FEVV


#include "FEVV/DataStructures/AIF/AIFFace.inl"
