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

#if defined _MSC_VER
#pragma warning(disable : 4251)
#endif

#include <iostream>
#include <vector>
#include <set>
#include <limits>
#include <stdexcept>
#include <boost/range/iterator_range.hpp>

#include "FEVV/DataStructures/AIF/AIFFace.hpp"
#include "FEVV/DataStructures/AIF/AIFVertex.hpp"
#include "FEVV/DataStructures/AIF/AIFEdge.hpp"

#include "FEVV/DataStructures/AIF/AIFProperties.h"
#include "FEVV/DataStructures/AIF/AIFCellContainer.h"

#include "FEVV/Tools/Container/Helpers.hxx"


namespace FEVV {
namespace DataStructures {
namespace AIF {

/**
 * \class	AIFMesh
 * \brief	This class represents an AIF structure.
 *			AIF structure can deal with both manifold and non-manifold
 *surface meshes. The faces of the meshes are not limited in degree, thus, the
 *structure also deals with polygonal meshes. To sum up, AIFMesh can handle
 *every surface mesh (and can be easily extended to upper dimensions). The
 *limitation of the structure is that the face with 0, 1 or 2 vertices are
 *unoriented (not complying case). For faces with at least three vertices,
 *orientation is deduced from its incident edge order. \see AIFVertex \see
 *AIFEdge \see		AIFFace
 */
class AIFMesh
{
public:
  typedef AIFMesh Self;
  typedef boost::shared_ptr< Self > ptr;
  typedef boost::shared_ptr< Self > ptr_mesh;
  typedef boost::shared_ptr< const Self > ptr_cmesh;

  typedef AIFVertex vertex_type;
  typedef AIFEdge edge_type;
  typedef AIFFace face_type;

  typedef AIFCellContainer< vertex_type::ptr > VertexContainerType;
  typedef AIFCellContainer< edge_type::ptr > EdgeContainerType;
  typedef AIFCellContainer< face_type::ptr > FaceContainerType;

  typedef AIFVertex::CoordinateType CoordinateType;
  typedef AIFPoint< CoordinateType, 3 > Point; // 3D point for geometry
  typedef std::array< CoordinateType, 2 >
      PointUV;                                   // 2D point for texture coords
  typedef AIFVector< CoordinateType, 3 > Vector; // 3D vector

  typedef AIFFace::NormalCoordinateType NormalCoordinateType;
  typedef AIFVector< NormalCoordinateType, 3 > Normal; // 3D vector

private:
  /*!
   *			Vertices (AIFVertex pointer) container of the mesh.
   */
  VertexContainerType m_Vertices;
  /*!
   *			Edges (AIFEdge pointer) container of the mesh.
   */
  EdgeContainerType m_Edges;
  /*!
   *			Faces (AIFFace pointer) container of the mesh.
   */
  FaceContainerType m_Faces;

  /*!
   *			Vertices property maps.
   */
  PropertyMapContainer m_VertexPropertyMaps;
  /*!
   *			Edges property maps.
   */
  PropertyMapContainer m_EdgePropertyMaps;
  /*!
   *			Faces property maps.
   */
  PropertyMapContainer m_FacePropertyMaps;
  /*!
   *			Associative property maps.
   */
  PropertyMapContainer m_AssocPropertyMaps;

public:
  /*!
   * 			Default constructor
   */
  AIFMesh(void);

  /*!
   * 			Copy constructor (deep copy)
   */
  AIFMesh(const Self &);

  /*!
  * 			Destructor
  */
  ~AIFMesh();

  /*!
   * 			= operator (deep copy)
   */
  Self &operator=(const Self &); 
public:
  /*!
   * 			Static constructor method
   * \return	A pointer to a new instance of AIFMesh class.
   */
  static ptr_mesh New();
  /*!
   * 			Mesh vertices getter
   * \return	An iterator range on the vertices (AIFVertex pointer) of the
   * involved mesh.
   */
  boost::iterator_range< VertexContainerType::const_iterator >
  GetVertices() const
  {
    return boost::make_iterator_range(m_Vertices.cbegin(), m_Vertices.cend());
  }
  /*!
   * 			Mesh vertices getter
   * \return	An iterator range on the vertices (AIFVertex pointer) of the
   * involved mesh.
   */
  boost::iterator_range< VertexContainerType::iterator > GetVertices()
  {
    return boost::make_iterator_range(m_Vertices.begin(), m_Vertices.end());
  }
  /*!
   * 			Mesh vertices begin iterator
   * \return	An iterator pointing on the beginning of the vertex container
   * (AIFVertex pointer) of the involved mesh.
   */
  inline VertexContainerType::iterator vertices_begin()
  {
    return m_Vertices.begin();
  }
  /*!
   * 			Mesh vertices end iterator
   * \return	An iterator pointing on the end of the vertex container
   * (AIFVertex pointer) of the involved mesh.
   */
  inline VertexContainerType::iterator vertices_end()
  {
    return m_Vertices.end();
  }

  /*!
   * 			Mesh edges getter
   * \return	An iterator range on the edges (AIFEdge pointer) of the involved
   * mesh.
   */
  boost::iterator_range< EdgeContainerType::const_iterator > GetEdges() const
  {
    return boost::make_iterator_range(m_Edges.cbegin(), m_Edges.cend());
  }
  /*!
   * 			Mesh edges getter
   * \return	An iterator range on the edges (AIFEdge pointer) of the involved
   * mesh.
   */
  boost::iterator_range< EdgeContainerType::iterator > GetEdges()
  {
    return boost::make_iterator_range(m_Edges.begin(), m_Edges.end());
  }
  /*!
   * 			Mesh faces getter
   * \return	An iterator range on the faces (AIFFace pointer) of the involved
   * mesh.
   */
  boost::iterator_range< FaceContainerType::const_iterator > GetFaces() const
  {
    return boost::make_iterator_range(m_Faces.cbegin(), m_Faces.cend());
  }
  /*!
   * 			Mesh faces getter
   * \return	An iterator range on the faces (AIFFace pointer) of the involved
   * mesh.
   */
  boost::iterator_range< FaceContainerType::iterator > GetFaces()
  {
    return boost::make_iterator_range(m_Faces.begin(), m_Faces.end());
  }
  /*!
   * 			Vertices count getter
   * \return	The number of vertices of the involved mesh.
   */
  unsigned int GetNumberOfVertices() const
  {
    return static_cast< unsigned int >(m_Vertices.size());
  }
  /*!
   * 			Edges count getter
   * \return	The number of edges of the involved mesh.
   */
  unsigned int GetNumberOfEdges() const
  {
    return static_cast< unsigned int >(m_Edges.size());
  }
  /*!
   * 			Faces count getter
   * \return	The number of faces of the involved mesh.
   */
  unsigned int GetNumberOfFaces() const
  {
    return static_cast< unsigned int >(m_Faces.size());
  }

  /*!
   * Insertion function for a simple vertex (without handling any
   * topological relation)
   *
   * \param  vertex The isolate vertex to add to the mesh
   */
  void InsertIsolatedVertex(vertex_type::ptr vertex) { m_Vertices.add(vertex); }
  /*!
   * Insertion function for a simple edge (without handling any
   * topological relation)
   *
   * \param  edge The isolate edge to add to the mesh
   */
  void InsertIsolatedEdge(edge_type::ptr edge) { m_Edges.add(edge); }
  /*!
   * Insertion function for a simple face (without handling any
   * topological relation)
   *
   * \param  face The isolate face to add to the mesh
   */
  void InsertIsolatedFace(face_type::ptr face) { m_Faces.add(face); }
  /*!
   * Suppression function for a simple cell (vertex, edge, face...) without
   * handling any topological relation.
   *
   * \param  cell The isolate cell to erase from the mesh
   */
  template< typename T >
  void EraseIsolatedCell(AIFCellContainer< T > &container, const T &cell)
  {
    std::size_t idx = cell->GetIndex();
    if (idx != -1)
    {
      // remove element from container
      std::size_t cLastId;
      if (container[idx]->GetIndex() == idx)
        cLastId = container.remove(idx);
      else
      {
        std::cerr << "EraseIsolatedCell: Warning: cell index does not match container index. A linear search is thus done." << std::endl;
        auto it = container.begin(), it_e = container.end();
        for (std::size_t cpt=0; it != it_e; ++it,++cpt)
        {
          if ((*it)->GetIndex() == idx)
          {
            cLastId = container.remove(cpt);
            break;
          }
        }
      }
      // remove entry from property maps
      PropertyMapContainer *pmc = GetPropertyMapContainer<T>();

      pmc->removeProperties(idx, cLastId);
      cell->SetIndex(-1);
    }
  }
  /*!
   * Suppression function for a simple vertex (without handling any
   * topological relation)
   *
   * \param  vertex The isolate vertex to erase from the mesh
   */
  void EraseIsolatedVertex(vertex_type::ptr vertex)
  {
    EraseIsolatedCell(m_Vertices, vertex);
  }
  /*!
   * Suppression function for a simple edge (without handling any
   * topological relation)
   * \param  edge The isolate edge to erase from the mesh
   */
  void EraseIsolatedEdge(edge_type::ptr edge)
  {
    EraseIsolatedCell(m_Edges, edge);
  }

  /*!
   * Suppression function for a simple edge iterator (without handling
   * any topological relation)
   *
   * \param  edge_iter The edge iterator to erase its associated isolated
   *                   edge from the mesh
   */
  void EraseIsolatedEdge(EdgeContainerType::iterator edge_iter)
  {
    m_Edges.erase(edge_iter);
  }
  /*!
   * Suppression function for a simple face (without handling any
   * topological relation)
   *
   * \param  face The isolate face to erase from the mesh
   */
  void EraseIsolatedFace(face_type::ptr face)
  {
    EraseIsolatedCell(m_Faces, face);
  }

  /*!
  * Removes all vertices, faces and edges from a graph.
  * Property maps associated with vertices, faces and edges 
  * are also cleared.
  */
  void clear();
public:
  /*!
   * Print the textual description of the involved mesh on the
   * standard output.
   */
  void Print() const;

  //---------------------------
  // Property maps management
  //---------------------------

  /*!
   *		  Get a property map container.
   *         This function is specialized for each kown CellType
   *         in AIFMesh.cpp.
   * \param  CellType  type of the cell the property map container
   *                   is associated to (AIFVertex/AIFEdge/AIFFace)
   */
  template< typename CellType >
  PropertyMapContainer *GetPropertyMapContainer(void);

  /*!
   *		  Get a property map container.
   *         This function is specialized for each kown CellType
   *         in AIFMesh.cpp.
   * \param  CellType  type of the cell the property map container
   *                   is associated to (AIFVertex/AIFEdge/AIFFace)
   */
  template< typename CellType >
  const PropertyMapContainer *GetPropertyMapContainer(void) const;

  /*!
   *		  Add a property map.
   * \param  CellType  type of the cell the property map container
   *                   is associated to (AIFVertex/AIFEdge/AIFFace)
   * \param  mapName  property map name
   */
  template< typename CellType, typename T >
  PropertyMap< T > *AddPropertyMap(const std::string &mapName)
  {
    PropertyMapContainer *pmc = GetPropertyMapContainer< CellType >();
    return pmc->addPropertyMap< T >(mapName);
  }

  /*!
   *		  Get a property map.
   * \param  CellType  type of the cell the property map container
   *                   is associated to (AIFVertex/AIFEdge/AIFFace)
   * \param  mapName  property map name
   */
  template< typename CellType, typename T >
  PropertyMap< T > *GetPropertyMap(const std::string &mapName)
  {
    PropertyMapContainer *pmc = GetPropertyMapContainer< CellType >();
    return pmc->getPropertyMap< T >(mapName);
  }


  /*!
   *		  Indicate if a property map exists for the given cell type.
   * \param  CellType  type of the cell the property map container
   *                   is associated to (AIFVertex/AIFEdge/AIFFace)
   * \param  mapName  property map name
   */
  template< typename CellType >
  bool isPropertyMap(const std::string &mapName) const
  {
    const PropertyMapContainer *pmc = GetPropertyMapContainer< CellType >();
    return pmc->isPropertyMap(mapName);
  }

  /*!
  *		  Indicate if at least one property map with given prefix exists for the given cell type.
  * \tparam  CellType  type of the cell the property map container
  *                   is associated to (AIFVertex/AIFEdge/AIFFace)
  * \param  mapPrefixName  property map prefix name
  */
  template< typename CellType >
  bool isAPropertyMapStartingWithPrefix(const std::string &mapPrefixName) const
  {
    const PropertyMapContainer *pmc = GetPropertyMapContainer< CellType >();
    return pmc->isAPropertyMapStartingWithPrefix(mapPrefixName);
  }

  /*!
  *		  Get property map names starting with prefix.
  * \tparam  CellType  type of the cell the property map container
  *                   is associated to (AIFVertex/AIFEdge/AIFFace)
  * \param  mapPrefixName  property map prefix name
  */
  template< typename CellType >
  std::vector< std::string > GetPropertyMapNamesStartingWithPrefix(const std::string &mapPrefixName) const
  {
    const PropertyMapContainer *pmc = GetPropertyMapContainer< CellType >();
    return pmc->GetPropertyMapNamesStartingWithPrefix(mapPrefixName);
  }

  /*!
   *		  Remove a property map.
   * \param  CellType  type of the cell the property map container
   *                   is associated to (AIFVertex/AIFEdge/AIFFace)
   * \param  mapName  property map name
   */
  template< typename CellType >
  void RemovePropertyMap(const std::string &mapName)
  {
    PropertyMapContainer *pmc = GetPropertyMapContainer< CellType >();
    pmc->removePropertyMap(mapName);
  }

  /*!
   *		  Set the property value for a given cell.
   * \param  CellType  type of the cell the property map container
   *                   is associated to (AIFVertex/AIFEdge/AIFFace)
   * \param  mapName  property map name
   * \param  cellId  cell (aka AIFVertex/AIFEdge/AIFFace) index
   * \param  value  value to set
   */
  template< typename CellType, typename T >
  void SetProperty(const std::string &mapName, std::size_t cellId, T value)
  {
    PropertyMapContainer *pmc = GetPropertyMapContainer< CellType >();
    pmc->setProperty< T >(mapName, cellId, value);
  }

  /*!
   *		  Get the property value of a given cell.
   * \param  CellType  type of the cell the property map container
   *                   is associated to (AIFVertex/AIFEdge/AIFFace)
   * \param  mapName  property map name
   * \param  cellId  cell (aka AIFVertex/AIFEdge/AIFFace) index
   */
  template< typename CellType, typename T >
  T &GetProperty(const std::string &mapName, std::size_t cellId)
  {
    PropertyMapContainer *pmc = GetPropertyMapContainer< CellType >();
    return pmc->getProperty< T >(mapName, cellId);
  }

  /*!
   *		  Get the property value of a given cell.
   * \param  CellType  type of the cell the property map container
   *                   is associated to (AIFVertex/AIFEdge/AIFFace)
   * \param  mapName  property map name
   * \param  cellId  cell (aka AIFVertex/AIFEdge/AIFFace) index
   */
  template< typename CellType, typename T >
  const T &GetProperty(const std::string &mapName, std::size_t cellId) const
  {
    const PropertyMapContainer *pmc = GetPropertyMapContainer< CellType >();
    return pmc->getProperty< T >(mapName, cellId);
  }

  //---------------------------
  // Associative property maps management
  //---------------------------

  /*!
   *         Add an associative property map.
   * \param  mapName    property map name
   */
  template< typename KeyType, typename ValueType >
  AssocPropertyMap< KeyType, ValueType > *
  AddAssocPropertyMap(const std::string &mapName)
  {
    return m_AssocPropertyMaps.addAssocPropertyMap< KeyType, ValueType >(
        mapName);
  }

  /*!
   *         Get an associative property map.
   * \param  mapName  property map name
   */
  template< typename KeyType, typename ValueType >
  AssocPropertyMap< KeyType, ValueType > *
  GetAssocPropertyMap(const std::string &mapName)
  {
    return m_AssocPropertyMaps.getAssocPropertyMap< KeyType, ValueType >(
        mapName);
  }

  /*!
   *         Remove an associative property map.
   * \param  mapName  property map name
   */
  void RemoveAssocPropertyMap(const std::string &mapName)
  {
    m_AssocPropertyMaps.removePropertyMap(mapName);
  }

  /*!
   *		  Indicate if an associative property map exists.
   * \param  mapName  property map name
   */
  bool isAssocPropertyMap(const std::string &mapName)
  {
    return m_AssocPropertyMaps.isPropertyMap(mapName);
  }

}; // END OF AIFMesh

template<>
PropertyMapContainer *
AIFMesh::GetPropertyMapContainer< AIFVertex::ptr >(void);
template<>
PropertyMapContainer *
AIFMesh::GetPropertyMapContainer< AIFEdge::ptr >(void);
template<>
PropertyMapContainer *
AIFMesh::GetPropertyMapContainer< AIFFace::ptr >(void);
template<>
const PropertyMapContainer *
AIFMesh::GetPropertyMapContainer< AIFVertex::ptr >(void) const;
template<>
const PropertyMapContainer *
AIFMesh::GetPropertyMapContainer< AIFEdge::ptr >(void) const;
template<>
const PropertyMapContainer *
AIFMesh::GetPropertyMapContainer< AIFFace::ptr >(void) const;

} // namespace AIF
} // namespace DataStructures
} // namespace FEVV


#include "FEVV/DataStructures/AIF/AIFMesh.inl"
