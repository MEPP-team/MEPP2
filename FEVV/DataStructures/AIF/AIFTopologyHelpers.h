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

#include <boost/range/iterator_range.hpp>
#include <exception>
#include <utility>
#include <algorithm>
#include <iostream>
#include <type_traits>

#include "FEVV/DataStructures/AIF/AIFMesh.hpp"

// The helpers functions are in 4 sections : vertex, edge, face, mesh. The
// functions are set in a section if the main parameter (thus the first)
// correspond to the section. eg : vertex_position(edge, vertex) for the
// position, V1 or V2, in an edge -> in the edge section.

// For the moment, only link_... and unlink_... functions realize topological
// changes (and it should stay like this until a really good reason) This
// statement certify that the topological caching (one ring, adjacency relations,
// ...) need to be delete only in these functions (the delete is actually realize
//by the private add_... and remove_... functions that are called)

// The topological caching is not yet refactored -> no one ring, adjacency
// relations, ... for now, at least with caching scheme

// For the moment, functions with no topological change use the static
// incident_...(vertex/edge/face) to get the incident relations whereas, the
// functions including topological changes use the non-const
// vertex/edge/face.Getincident...() functions

namespace FEVV {
namespace DataStructures {
namespace AIF {

/**
 * \class	AIFTopologyHelpers
 * \brief	This class is an helper class associated to the AIFMesh
 *structure. AIFTopologyHelpers implements all the basic topologic function used
 *to manipulate an AIFMesh. Any operation on the AIF structure need to be
 *realized with this helper class. \see		AIFMesh
 */
class AIFTopologyHelpers
{
public:
  typedef AIFMesh mesh_type;
  typedef mesh_type::vertex_type vertex_type;
  typedef mesh_type::edge_type edge_type;
  typedef mesh_type::face_type face_type;

  typedef mesh_type::ptr smart_ptr_mesh;
  typedef mesh_type::ptr_cmesh smart_ptr_cmesh;
  typedef mesh_type *ptr_mesh;
  typedef const mesh_type *ptr_cmesh;
  typedef mesh_type &ref_mesh;
  typedef typename std::add_lvalue_reference< const mesh_type >::type ref_cmesh;
  typedef vertex_type::ptr vertex_descriptor;
  typedef edge_type::ptr edge_descriptor;
  typedef face_type::ptr face_descriptor;

  typedef vertex_type::VertexContainerType vertex_container_in_vertex;
  typedef vertex_type::EdgeContainerType edge_container_in_vertex;
  typedef vertex_type::FaceContainerType face_container_in_vertex;
  typedef edge_type::VertexContainerType vertex_container_in_edge;
  typedef edge_type::EdgeContainerType edge_container_in_edge;
  typedef edge_type::FaceContainerType face_container_in_edge;
  typedef face_type::VertexContainerType vertex_container_in_face;
  typedef face_type::EdgeContainerType edge_container_in_face;
  typedef face_type::FaceContainerType face_container_in_face;
  typedef mesh_type::VertexContainerType vertex_container;
  typedef mesh_type::EdgeContainerType edge_container;
  typedef mesh_type::FaceContainerType face_container;

  typedef unsigned int size_type;

  typedef edge_container_in_vertex::const_iterator out_edge_iterator;
  typedef edge_container_in_vertex::size_type degree_size_type;

  typedef vertex_container_in_vertex::const_iterator adjacency_iterator;

  typedef edge_container_in_face::iterator incident_edge_iterator;
  typedef face_container_in_edge::iterator incident_face_iterator;

  typedef edge_container_in_face::const_iterator const_incident_edge_iterator;
  typedef face_container_in_edge::const_iterator const_incident_face_iterator;

  /*!
   * \enum		vertex_pos
   * \brief	The vertex position in an edge
   */
  enum vertex_pos { FIRST, SECOND };
  /*!
   * \enum		face_degree
   * \brief	Compliant enumeration for face degree comparison
   */
  enum face_degree { TRIANGULAR = 3, QUADRANGULAR = 4 };

  // nulls
  static vertex_descriptor null_vertex() { return vertex_descriptor(); }
  static edge_descriptor null_edge() { return edge_descriptor(); }
  static face_descriptor null_face() { return face_descriptor(); }

public:
  //////////////////////////////////	Vertex AIFTopologyHelpers
  /////////////////////////////////////
  /*!
   * Function determining if the argument vertex is degenerated
   *
   * \param  vertex  The involving vertex
   * \return return  true if vertex==null_vertex, and false otherwise.
   */
  static bool is_degenerated_vertex(vertex_descriptor vertex)
  {
    return vertex == null_vertex();
  }

  /*!
   * 			Vertex degree
   * \param	vertex	The involving vertex
   * \return	The number of incident edges of the argument vertex.
   */
  static size_type degree(vertex_descriptor vertex)
  {
    return vertex->GetDegree();
  }
  /*!
   * Function determining if the argument vertex is isolated (without
   * any incident edge)
   *
   * \param  vertex  The involving vertex
   * \return true if the argument vertex is isolated, false otherwise.
   */
  static bool is_isolated_vertex(vertex_descriptor vertex)
  {
    return vertex->GetDegree() == 0u;
  }
  /*!
   * Function determining if the argument vertex is a cut vertex
   *
   * \param  vertex  The involving vertex
   * \return true if the argument vertex is a cut vertex, false otherwise.
   */
  static bool is_cut_vertex(vertex_descriptor vertex)
  {
    typedef edge_container_in_vertex::const_iterator it_type;
    int nbBoundaryEdges = 0;
    boost::iterator_range< it_type > edges_range = incident_edges(vertex);

    for(it_type it = edges_range.begin(); it != edges_range.end(); ++it)
    {
      if(is_surface_border_edge(*it))
        nbBoundaryEdges++;
    }

    if((nbBoundaryEdges < 4) || (nbBoundaryEdges % 2 == 1))
      return false;
    else
      return ((nbBoundaryEdges - 4) % 2 == 0);
  }

  /*!
   * Function determining if the argument vertex is a 2-manifold vertex
   * considering only the topology (not the geometry)
   *
   * \param  vertex  The involving vertex
   *
   * \return  true if the argument vertex is a locally 2-manifold vertex,
              false otherwise.
   */
  static bool is_2_manifold_vertex(vertex_descriptor vertex)
  {
    if(is_degenerated_vertex(vertex) || is_isolated_vertex(vertex) ||
       is_cut_vertex(vertex))
      return false;

    auto incidentEdgesRange = incident_edges(vertex);
    auto it = incidentEdgesRange.begin();
    auto ite = incidentEdgesRange.end();
    for(; it != ite; ++it)
    {
      if(is_degenerated_edge(*it) || is_isolated_edge(*it) ||
         is_dangling_edge(*it) || is_complex_edge(*it))
        return false;
    }

    return true;
  }

  /*!
   * Function determining if the argument vertex is on a border (one of
   * its incident edges is on border)
   *
   * \param  vertex  The involving vertex
   *
   * \return  true if the argument vertex is on a border, false otherwise.
   */
  static bool is_surface_border_vertex(vertex_descriptor vertex)
  {
    if(is_isolated_vertex(vertex))
      return false;

    typedef edge_container_in_vertex::const_iterator it_type;
    boost::iterator_range< it_type > edges_range = incident_edges(vertex);

    for(it_type it = edges_range.begin(); it != edges_range.end(); ++it)
    {
      if(is_surface_border_edge(*it))
        return true;
    }
    return false;
  }
  /*!
   * Function determining if the argument vertex is not on a border
   * (none of its incident edges is on border)
   *
   * \param  vertex  The involving vertex
   *
   * \return  true if the argument vertex is not on a border and not an
   *          isolated vertex, false otherwise.
   */
  static bool is_surface_interior_vertex(vertex_descriptor vertex)
  {
    return !is_surface_border_vertex(vertex) && !is_isolated_vertex(vertex);
  }
  /*!
   * \brief   Function determining if the arguments share an incidence relation
   *
   * \param   vertex  The involving vertex
   * \param   edge    The involving edge
   *
   * \return  true if the arguments share an incidence relation,
   *          false otherwise.
   * \warning are_incident(edge, vertex) will not give consistent
   *          results when topological modifications are in progress
   */
  static bool are_incident(vertex_descriptor vertex, edge_descriptor edge)
  {
	if(vertex==null_vertex())
		return false;	  
    typedef edge_container_in_vertex::const_iterator it_type;
    boost::iterator_range< it_type > edges_range = incident_edges(vertex);

    for(it_type it = edges_range.begin(); it != edges_range.end(); ++it)
    {
      if(*it == edge)
        return true;
    }
    return false;
  }
  /*!
   * \brief  Function determining if the arguments share an incidence relation
   *
   * \param  vertex  The involving vertex
   * \param  u       one vertex of the involving edge
   * \param  v       the other vertex of the involving edge
   *
   * \return true if the arguments share an incidence relation, false
   *         otherwise.
   */
  static bool are_incident(vertex_descriptor vertex,
                           vertex_descriptor u,
                           vertex_descriptor v)
  {
	if(vertex==null_vertex())
		return false;
	  
    edge_descriptor edge = common_edge(u, v);
    if(edge != null_edge())
    {
      typedef edge_container_in_vertex::const_iterator it_type;
      boost::iterator_range< it_type > edges_range = incident_edges(vertex);

      for(it_type it = edges_range.begin(); it != edges_range.end(); ++it)
      {
        if(*it == edge)
          return true;
      }
    }
    return false;
  }
  /*!
   * Function determining if the arguments share an adjacent relation
   * (if they share an edge)
   *
   * \param  vertex1  The first involving vertex
   * \param  vertex2  The second involving vertex
   *
   * \return  true if the arguments share an adjacency relation, false
   *          otherwise (a vertex is not adjacent to itself).
   */
  static bool are_adjacent(vertex_descriptor vertex1, vertex_descriptor vertex2)
  {
    typedef edge_container_in_vertex::const_iterator it_type;
    boost::iterator_range< it_type > edges_range = incident_edges(vertex1);

    for(it_type it = edges_range.begin(); it != edges_range.end(); ++it)
    {
      if(opposite_vertex(*it, vertex1) == vertex2)
        return true;
    }
    return false;
  }
  /*!
   * 			Adjacent relations with vertices
   * \param	vertex	The involving vertex
   * \return	An iterator range on the adjacent vertices (AIFVertex pointer) of
   * the involved vertex.
   */
  static std::vector< vertex_descriptor >
  adjacent_vertices(vertex_descriptor vertex)
  {
    return get_ordered_one_ring_of_adjacent_vertices(vertex);
  }
  /*!
   * 			Incidence relations with edges
   * \param	vertex	The involving vertex
   * \return	An iterator range on the incident edges (AIFEdge pointer) of the
   * involved vertex.
   */
  static boost::iterator_range< edge_container_in_vertex::const_iterator >
  incident_edges(vertex_descriptor vertex)
  {
    if(vertex == null_vertex())
      throw std::invalid_argument(
          "AIFTopologyHelpers::incident_edges -> no incident relation exists "
          "for a null vertex.");

    return vertex->GetIncidentEdges();
  }
  /*!
   * 			Incidence relations with faces
   * \param	vertex	The involving vertex
   * \return	An iterator range on the incident faces (AIFFace pointer) of the
   * involved vertex.
   */
  static boost::iterator_range< face_container_in_vertex::const_iterator >
  incident_faces(vertex_descriptor vertex)
  {
    if(!vertex->m_Incident_PtrFaces_Computed)
    {
      typedef edge_container_in_vertex::const_iterator it_type;
      typedef face_container_in_edge::const_iterator it_type2;

      boost::iterator_range< it_type > edges_range = incident_edges(vertex);

      std::set< face_descriptor > incFaces;

      for(it_type it = edges_range.begin(); it != edges_range.end(); ++it)
      {
        boost::iterator_range< it_type2 > faces_range =
            incident_faces(*it); // get all incident faces of current edge
        incFaces.insert(faces_range.begin(), faces_range.end());
      }
      face_container_in_vertex common_faces(incFaces.begin(), incFaces.end());
      vertex->m_Incident_PtrFaces = common_faces;
      vertex->m_Incident_PtrFaces_Computed = true;
    }
    return boost::make_iterator_range(vertex->m_Incident_PtrFaces.cbegin(),
                                      vertex->m_Incident_PtrFaces.cend());
  }
  /*!
   * 			Incidence relations with faces
   * \param	vertex	The involving vertex
   * \return	A container of the incident faces (AIFFace pointer) of the
   * involved vertex.
   */
  static face_container_in_vertex
  incident_faces_container(vertex_descriptor vertex)
  {
    incident_faces(vertex); // update calculus
    return face_container_in_vertex(vertex->m_Incident_PtrFaces.cbegin(),
                                    vertex->m_Incident_PtrFaces.cend());
  }
  /*!
   * 			First common edge between two vertices
   * \param	vertex1	The first involving vertex
   * \param	vertex2	The second involving vertex
   * \return	The first common edge between vertex1 and vertex2.
   */
  static edge_descriptor common_edge(vertex_descriptor vertex1,
                                     vertex_descriptor vertex2)
  {
    typedef edge_container_in_vertex::const_iterator it_type;
    // if (vertex1 == null_vertex() || vertex2 == null_vertex())
    //  return null_edge();
    std::vector< edge_descriptor > common_edges;
    boost::iterator_range< it_type > edges_range = incident_edges(
        vertex1); // we need to use the first vertex (the source in concept!!)


    for(it_type it = edges_range.begin(); it != edges_range.end(); ++it)
    {
      if(opposite_vertex(*it, vertex1) == vertex2)
        return *it;
    }
    return null_edge();
  }
  /*!
   * 			Common edges between two vertices
   * \param	vertex1	The first involving vertex
   * \param	vertex2	The second involving vertex
   * \return	An std::vector of the common edges between vertex1 and vertex2
   * (there can be 0, 1, or n common edges).
   */
  static std::vector< edge_descriptor > common_edges(vertex_descriptor vertex1,
                                                     vertex_descriptor vertex2)
  {
    typedef edge_container_in_vertex::const_iterator it_type;
    std::vector< edge_descriptor > common_edges;
    boost::iterator_range< it_type > edges_range = incident_edges(vertex1);

    for(it_type it = edges_range.begin(); it != edges_range.end(); ++it)
    {
      if(opposite_vertex(*it, vertex1) == vertex2)
        common_edges.push_back(*it);
    }
    return common_edges;
  }
  /*!
   * Vertex and edge link.
   * An incidence relation is created between the vertex and the edge
   * by adding the vertex at the specified position in the edge, and adding the
   * edge in the incident edges container of the vertex. Update caching
   * information for vertices.
   *
   * \param  vertex    The involving vertex
   * \param  edge      The involving edge
   * \param  position  The position of the vertex in the edge (FIRST or
   *                   SECOND)
   *
   * \warning  Be careful, this function will replace the former vertex
   *           at the position \c position in the edge. Thus you HAVE to
   *           unlink these two element before, otherwise topological
   *           incoherences will be introduce.
   * \see      add_vertex_to_edge(edge, vertex, position)
   * \see      add_edge_to_vertex(vertex, edge)
   */
  static void link_vertex_and_edge(vertex_descriptor vertex,
                                   edge_descriptor edge,
                                   vertex_pos position)
  {
    if(!are_incident(vertex, edge))
    {
      add_vertex_to_edge(edge, vertex, position);
      add_edge_to_vertex(vertex, edge);
    }
    else
      throw std::invalid_argument(
          "AIFTopologyHelpers::link_vertex_and_edge -> incident relation "
          "already exists between the two elements.");
  }

  /*!
   * Vertex and edge unlink.
   * An incidence relation is removed between the vertex and the edge
   * by removing the vertex at the specified position in the edge, and removing
   * the edge in the incident edges container of the vertex.
   *
   * \param  vertex  The involving vertex
   * \param  edge    The involving edge
   *
   * \warning  Be careful, this function leave the edge with only one vertex
   *           (the other become a null_vertex)
   * \see      remove_vertex_from_edge(edge, vertex)
   * \see      remove_edge_from_vertex(vertex, edge)
   */
  static void unlink_vertex_and_edge(vertex_descriptor vertex,
                                     edge_descriptor edge)
  {
    bool edge_has_that_incident_vertex = are_incident(edge, vertex);
    bool vertex_has_that_incident_edge = are_incident(vertex, edge);

    if(!edge_has_that_incident_vertex && !vertex_has_that_incident_edge)
      throw std::invalid_argument(
          "AIFTopologyHelpers::unlink_vertex_and_edge -> no incident relation "
          "existing between the two elements.");

    if(edge_has_that_incident_vertex)
    {
      remove_vertex_from_edge(edge, vertex);
    }

    if(vertex_has_that_incident_edge)
    {
      remove_edge_from_vertex(vertex, edge);
    }
  }

  /*!
  * 			Sort mesh vertices
  * \param	mesh	The involving mesh
  * \param	cmp   The involving comparator object used for the sorting. Comparator can be based
  *               on natural ordering, on a spanning tree, etc.
  * \return	An iterator range on the vertices (AIFVertex pointer) of the involving mesh.
  */
  template<typename ComparatorType>
  static boost::iterator_range<vertex_container::const_iterator> sort_vertices(ptr_mesh mesh, ComparatorType cmp)
  {
    std::sort(mesh->GetVertices().begin(), mesh->GetVertices().end(), cmp);
    return mesh->GetVertices();
  }

  /*!
  * 			Sort mesh vertices
  * \param	mesh	The involving mesh
  * \param	cmp   The involving comparator object used for the sorting. Comparator can be based
  *               on natural ordering, on a spanning tree, etc.
  * \return	An iterator range on the vertices (AIFVertex pointer) of the involving mesh.
  */
  template<typename ComparatorType>
  static boost::iterator_range<vertex_container::const_iterator> sort_vertices(smart_ptr_mesh mesh, ComparatorType cmp)
  {
    return sort_vertices<ComparatorType>(mesh.get(), cmp);
  }

  /*!
  * 			Sort mesh vertices
  * \param	mesh	The involving mesh
  * \param	cmp   The involving comparator object used for the sorting. Comparator can be based
  *               on natural ordering, on a spanning tree, etc.
  * \return	An iterator range on the vertices (AIFVertex pointer) of the involving mesh.
  */
  template<typename ComparatorType>
  static boost::iterator_range<vertex_container::const_iterator> sort_vertices(ref_mesh mesh, ComparatorType cmp)
  {
    return sort_vertices<ComparatorType>(&mesh, cmp);
  }
  /*!
  * 			Sort incident edges starting with given incident edge and return incidence relations with edges
  * \param	vertex	The involving vertex
  * \param	cmp   	The involving comparator object used for the sorting. Comparator can be based
  *                 on natural ordering, on a spanning tree, on one-ring vertex order (if incident
  *                 faces are counterclokwise oriented, then one-ring is clockwise oriented) etc.
  * \param	edge  	The involving incident edge to vertex
  * \param  do_full_incident_edge_sorting A boolean to do a sorting before starting at edge edge.
                    Sometimes we do no want to lose the current incident ordering (e.g. clockwise).
					In this case, we need to specify false.
  * \return	An iterator range on the incident edges (AIFEdge pointer) of the involve vertex.
  */
  template<typename ComparatorType>
  static boost::iterator_range<edge_container_in_vertex::const_iterator> sort_incident_edges_starting_with_edge(vertex_descriptor vertex, 
                                                                                                                ComparatorType cmp, 
                                                                                                                edge_descriptor edge,
                                                                                                                bool do_full_incident_edge_sorting)
  {
    if(do_full_incident_edge_sorting)
      std::sort(vertex->m_Incident_PtrEdges.begin(), vertex->m_Incident_PtrEdges.end(), cmp); // possible break of clockwise ordering (but no matter since usual order is insertion order)
    /////////////////////////////////////////////////////////////////////////
    std::list<edge_descriptor> ltmp, ltmpnext;
    edge_container_in_vertex::const_iterator it = vertex->m_Incident_PtrEdges.begin(),
      ite = vertex->m_Incident_PtrEdges.end();
    while (it != ite && *it != edge)
    {
      ltmpnext.push_back(*it);
      ++it;
    }
    while (it != ite)
    {
      ltmp.push_back(*it);
      ++it;
    }
    /////////////////////////////////////////////////////////////////////////
    vertex->m_Incident_PtrEdges.clear();
    vertex->m_Incident_PtrEdges.insert(vertex->m_Incident_PtrEdges.end(), ltmp.begin(), ltmp.end());
    vertex->m_Incident_PtrEdges.insert(vertex->m_Incident_PtrEdges.end(), ltmpnext.begin(), ltmpnext.end());
    /////////////////////////////////////////////////////////////////////////
    return incident_edges(vertex);
  }

  /*!
  * 			Sort incident edges and return incidence relations with edges
  * \param	vertex	The involving vertex
  * \param	cmp   	The involving comparator object used for the sorting. Comparator can be based
  *                 on natural ordering, on a spanning tree, on one-ring vertex order (if incident
  *                 faces are counterclokwise oriented, then one-ring is clockwise oriented) etc.
  * \param  do_full_incident_edge_sorting A boolean to do a sorting before starting at minimun edge.
  * \return	An iterator range on the incident edges (AIFEdge pointer) of the involve vertex.
  */
  template<typename ComparatorType>
  static boost::iterator_range<edge_container_in_vertex::const_iterator> sort_incident_edges( vertex_descriptor vertex, 
                                                                                              ComparatorType cmp,
                                                                                              bool do_full_incident_edge_sorting)
  {

    return sort_incident_edges_starting_with_edge<ComparatorType>(vertex, 
                                                                  cmp, 
                                                                  *std::min_element(vertex->m_Incident_PtrEdges.begin(), 
                                                                                    vertex->m_Incident_PtrEdges.end(), 
                                                                                    cmp),
                                                                  do_full_incident_edge_sorting);
  }

  /*!
  * 			Sort incident edges of its two vertices
  * \param	edge	The involving edge
  * \param	cmp   The involving comparator object used for the sorting. Comparator can be based
  *               on natural ordering, on a spanning tree, on one-ring vertex order (if incident
  *               faces are counterclokwise oriented, then one-ring is clockwise oriented) etc.
  * \param  do_full_incident_edge_sorting A boolean to do a sorting before starting at minimun edge.
  */
  template<typename ComparatorType>
  static void sort_incident_edges(edge_descriptor edge,
                                  ComparatorType cmp,
                                  bool do_full_incident_edge_sorting)
  {
    sort_incident_edges<ComparatorType>(edge->get_first_vertex(), cmp, do_full_incident_edge_sorting);
    sort_incident_edges<ComparatorType>(edge->get_second_vertex(), cmp, do_full_incident_edge_sorting);
  }

  /*!
  * 			Sort incident edges of its two vertices' one-ring
  * \param	edge	The involving edge
  * \param	cmp   The involving comparator object used for the sorting. Comparator can be based
  *               on natural ordering, on a spanning tree, on one-ring vertex order (if incident
  *               faces are counterclokwise oriented, then one-ring is clockwise oriented) etc.
  */
  template<typename ComparatorType>
  static void sort_one_ring_incident_edges( edge_descriptor edge,
                                            ComparatorType cmp)
  {
    typedef edge_container_in_vertex::const_iterator it_type;
    boost::iterator_range<it_type> edges_range = incident_edges(edge->get_first_vertex());
    for (it_type it = edges_range.begin(); it != edges_range.end(); ++it) 
      sort_incident_edges<ComparatorType>(opposite_vertex(*it, edge->get_first_vertex()), cmp);

    edges_range = incident_edges(edge->get_second_vertex());
    for (it_type it = edges_range.begin(); it != edges_range.end(); ++it)
      sort_incident_edges<ComparatorType>(opposite_vertex(*it, edge->get_second_vertex()), cmp);

  }

  
  /*!
  *      Get incident border edges within a mesh hole. 
  * \param	 vertex	The involving vertex
  * \return  An std::vector of incident hole border edges. 
  * \warning vertex is assumed to be on the hole border, and no check
  *          is done.
  */
  static std::vector<edge_descriptor> get_incident_hole_border_edges(vertex_descriptor vertex)
  {
    std::vector<edge_descriptor> res;
    if (is_isolated_vertex(vertex)) return res;

    typedef edge_container_in_vertex::const_iterator it_type;
    boost::iterator_range<it_type> edges_range = incident_edges(vertex);

    for (it_type it = edges_range.begin(); it != edges_range.end(); ++it) {
      if ((is_surface_border_edge(*it) && 
          (degree((*it)->get_first_vertex())>2) && 
          (degree((*it)->get_second_vertex())>2)) ||
          (is_dangling_edge(*it) && 
          (degree((*it)->get_first_vertex()) >= 2) && 
          (degree((*it)->get_second_vertex()) >= 2))
        )
        res.push_back(*it);
    }
    if (res.size() == 1)
    {
      std::vector<edge_descriptor> res_init(res.begin(), res.end());
      for (it_type it = edges_range.begin(); it != edges_range.end(); ++it) {
        if (std::find(res_init.begin(), res_init.end(), *it) == res_init.end())
        {
          auto iter = res_init.begin();
          for (; iter != res_init.end(); ++iter)
          {
            if (are_adjacent(opposite_vertex(*it, vertex), opposite_vertex(*iter, vertex)))
            {
              res.push_back(*it);
              break;
            }
          }
        }
      }
    }
    return res;
  }

  /*!
  * 			Return a vector of incident edges to vertex except edge_to_remove
  * \param	vertex	The involving vertex
  * \param	edge_to_remove    The involving edge to remove from incident edges
  * \return	A std::vector of incident edges to vertex except edge_to_remove.
  * \see get_incident_hole_border_edges
  */
  static std::vector<edge_descriptor> get_incident_hole_border_edges_except_one_edge(vertex_descriptor vertex, edge_descriptor edge_to_remove)
  {
    std::vector<edge_descriptor> res = get_incident_hole_border_edges(vertex);

    auto iter = std::find(res.begin(), res.end(), edge_to_remove);
    if (iter != res.end())
      res.erase(iter);

    return res;
  }

  /*!
  * 			First common edge between a vertex and a face
  * \param	vertex	The involving vertex
  * \param	face    The involving face
  * \return	The first common edge between vertex and face.
  */
  static  edge_descriptor common_edge(vertex_descriptor vertex, face_descriptor face)
  {
    typedef edge_container_in_vertex::const_iterator it_type;

    boost::iterator_range<it_type> edges_range = incident_edges(vertex); // we need to use the first vertex (the source in concept!!)

    for (it_type it = edges_range.begin(); it != edges_range.end(); ++it) {
      if (are_incident(face, *it))
        return *it;
    }
    return null_edge();
  }

  /*!
  * 			First common edge between a vertex and a face except not_that_edge
  * \param	vertex	The involving vertex
  * \param	face    The involving face
  * \param	not_that_edge The involving edge
  * \return	The first common edge between vertex and face except not_that_edge.
  */
  static  edge_descriptor common_edge(vertex_descriptor vertex, face_descriptor face, edge_descriptor not_that_edge)
  {
    typedef edge_container_in_vertex::const_iterator it_type;

    boost::iterator_range<it_type> edges_range = incident_edges(vertex); // we need to use the first vertex (the source in concept!!)

    for (it_type it = edges_range.begin(); it != edges_range.end(); ++it) {
      if (are_incident(face, *it) && (*it != not_that_edge))
        return *it;
    }
    return null_edge();
  }
  //////////////////////////////////	Edge AIFTopologyHelpers
  /////////////////////////////////////
  /*!
   * 			Edge degree
   * \param	edge	The involving edge
   * \return	The number of incident faces of the argument edge.
   */
  static size_type degree(edge_descriptor edge) { return edge->GetDegree(); }
  /*!
   * 			Function determining if the argument edge is on a border
   * \param	edge	The involving edge
   * \return	true if the argument edge has only one incident face, false
   * otherwise. Note that dangling edges are NOT border edges.
   */
  static bool is_surface_border_edge(edge_descriptor edge)
  {
    return degree(edge) == 1;
  }
  /*!
   * Function determining if the argument edge is not on a border
   *
   * \param  edge  The involving edge
   *
   * \return  true if the argument edge is not on a border, false otherwise.
   */
  static bool is_surface_interior_edge(edge_descriptor edge)
  {
    return !is_surface_border_edge(edge) && !is_dangling_edge(edge);
  }

  /*!
   * 			Function determining if the argument edge is isolated
   *       (without any incident face + its vertices are 1 degree vertices)
   * \param	edge	The involving edge
   * \return	true if the argument edge is isolated, false otherwise.
   */
  static bool is_isolated_edge(edge_descriptor edge)
  {
    return (edge->get_first_vertex()->GetDegree() == 1u) &&
           (edge->get_second_vertex()->GetDegree() == 1u);
  }

  /*!
   * Function determining if the argument edge is a dangling edge
   *
   * \param  edge  The involving edge
   *
   * \return  true if the argument edge does not have any incident face and
   *          is not isolated, false otherwise.
   */
  static bool is_dangling_edge(edge_descriptor edge)
  {
    return (edge->GetDegree() == 0u) && !is_isolated_edge(edge);
  }
  /*!
   * Function determining if the argument edge is a complex edge
   *
   * \param  edge  The involving edge
   *
   * \return  true if the argument edge has strictly more than two incident
   *          face, false otherwise.
   */
  static bool is_complex_edge(edge_descriptor edge)
  {
    return edge->GetDegree() > 2u;
  }
  /*!
   * 			Function determining if the argument edge is degenerated
   * \param	edge	The involving edge
   * \return	Noting v1 and v2 the two incident vertices of the argument edge,
   *the function return true if edge is null or v1==null_vertex or
   *v2==null_vertex or v1==v2, and false otherwise. \warning  v1==v2 may occur
   *during a topological modification in progress
   */
  static bool is_degenerated_edge(edge_descriptor edge)
  {
    return edge == null_edge() || // a null edge is degenerated
           edge->get_first_vertex() ==
               edge->get_second_vertex() || // null-length edge are
                                            // degenereated, thus edge with 2
                                            // identical vertices are too.
           is_degenerated_vertex(edge->get_first_vertex()) ||
           is_degenerated_vertex(edge->get_second_vertex());
  }

  /*!
   * Function determining the common vertex between 2 incident edges.
   * Usage precondition: edge1 and edge2 must not be similar/parallel
   * (sharing the 2 same vertices).
   *
   * \param  edge1  The first involved edge
   * \param  edge2  The second involved edge
   *
   * \return  A vertex descriptor, this descriptor being a null vertex if the 2
   *         edges do not share any vertex.
   */
  static vertex_descriptor common_vertex(edge_descriptor edge1,
                                         edge_descriptor edge2)
  {
    if(edge1 == edge2)
      throw std::invalid_argument("AIFTopologyHelpers::common_vertex -> the 2 "
                                  "edge arguments are identical.");

    if((edge1->get_first_vertex() == edge2->get_first_vertex()) ||
       (edge1->get_first_vertex() == edge2->get_second_vertex()))
    {
      return edge1->get_first_vertex();
    }
    else if((edge1->get_second_vertex() == edge2->get_first_vertex()) ||
            (edge1->get_second_vertex() == edge2->get_second_vertex()))
      return edge1->get_second_vertex();

    return null_vertex();
  }
  
	/*!
	* 			Function counting the number of incident dangling edges
	* \param	vertex	The involving vertex
	* \return	number of incident dangling edges.
	*/
	static size_type num_incident_dangling_edge(vertex_descriptor vertex)
	{
	  size_type cpt = 0;
	  auto edges_range = incident_edges(vertex);
	  auto iterE = edges_range.begin();
	  for (; iterE != edges_range.end(); ++iterE)
	  {
		if (is_dangling_edge(*iterE))
		  ++cpt;
	  }
	  return cpt;
	}
	/*!
	* 			Function determining if the argument vertex has at least one incident dangling edge
	* \param	vertex	The involving vertex
	* \return	true if the argument vertex has at least one incident dangling edge, false otherwise.
	*/
	static bool has_incident_dangling_edge(vertex_descriptor vertex)
	{
	  return (num_incident_dangling_edge(vertex)>0);
	}

	/*!
	* 			Function returning the first incident dangling edge
	* \param	vertex	The involving vertex
	* \return	the first incident dangling edge if any, else return null_edge descriptor.
	*/
	static edge_descriptor get_incident_dangling_edge(vertex_descriptor vertex)
	{
	  edge_descriptor de = null_edge();
	  auto edges_range = incident_edges(vertex);
	  auto iterE = edges_range.begin();
	  for (; iterE != edges_range.end(); ++iterE)
	  {
		if (is_dangling_edge(*iterE))
		{
		  de = *iterE;
		  break;
		}
	  }
	  return de;
	}  
  /*!
   * \brief  Function determining if the argument edge shares an incidence
   *         relation with the argument vertex.
   *
   * \param  edge    The involving edge
   * \param  vertex  The involving vertex
   *
   * \return true if the arguments share an incidence relation, false
   *         otherwise.
   */
  static bool are_incident(edge_descriptor edge, vertex_descriptor vertex)
  {
	if(edge == null_edge())
      return false;		
    return (edge->get_first_vertex() == vertex) ||
           (edge->get_second_vertex() == vertex);
  }
  /*!
   * \brief  Function determining if the argument edge share an incidence
   *         relation with the argument face.
   *
   * \param  edge  The involving edge
   * \param  face  The involving face
   *
   * \return true if the arguments share an incidence relation, false
   *         otherwise.
   */
  static bool are_incident(edge_descriptor edge, face_descriptor face)
  {
	if(edge == null_edge())
      return false;		  
    typedef face_container_in_edge::const_iterator it_type;
    boost::iterator_range< it_type > faces_range = incident_faces(edge);

    for(it_type it = faces_range.begin(); it != faces_range.end(); ++it)
    {
      if(*it == face)
        return true;
    }
    return false;
  }
  /*!
   * 			Function determining if the arguments share an adjacent relation
   * (if they share a vertex) Usage precondition: edge1 and edge2 must not be
   * similar/parallel (sharing the 2 same vertices).
   * \param	edge1	The first involving edge
   * \param	edge2	The second involving edge
   * \return	true if the arguments share an adjacent relation, false otherwise
   * (an edge is not adjacent to itself).
   */
  static bool are_adjacent(edge_descriptor edge1, edge_descriptor edge2)
  {
    if(edge1 == edge2)
      return false;
    return are_incident(edge2, edge1->get_first_vertex()) ||
           are_incident(edge2, edge1->get_second_vertex());
  }

  /*!
   * 			Incidence relations with vertices
   * \param	edge	The involving edge
   * \return	An iterator range on the incident vertices (AIFVertex pointer) of
   * the involved edge.
   */
  static boost::iterator_range< vertex_container_in_edge::const_iterator >
  incident_vertices(edge_descriptor edge)
  { // Return type can't be iterator_range type cause vertex_container_in_edge
    // is an std::pair
    if(edge == null_edge())
      throw std::invalid_argument("AIFTopologyHelpers::incident_vertices -> no "
                                  "incident relation exists for a null edge.");
    return edge->GetIncidentVertices();
  }
  /*!
   * 			Adjacent relations with edges
   * \param	edge	The involving edge
   * \return	An iterator range on the adjacent edges (AIFEdge pointer) of the
   * involved edge.
   */
  static std::vector< edge_descriptor > adjacent_edges(edge_descriptor edge)
  {
    std::vector< edge_descriptor > adEdges1(
        edge->get_first_vertex()->GetIncidentEdges().begin(),
        edge->get_first_vertex()->GetIncidentEdges().end());
    std::vector< edge_descriptor > adEdges2(
        edge->get_second_vertex()->GetIncidentEdges().begin(),
        edge->get_second_vertex()->GetIncidentEdges().end());

    std::vector< edge_descriptor > res(edge->get_first_vertex()->GetDegree() +
                                       edge->get_second_vertex()->GetDegree() -
                                       2);

    std::sort(adEdges1.begin(), adEdges1.end());
    std::sort(adEdges2.begin(), adEdges2.end());
    adEdges1.erase(std::find(adEdges1.begin(), adEdges1.end(), edge));
    adEdges2.erase(std::find(adEdges2.begin(), adEdges2.end(), edge));

    std::set_union(adEdges1.begin(),
                   adEdges1.end(),
                   adEdges2.begin(),
                   adEdges2.end(),
                   res.begin());

    return res; // note that it is a non-sens to manage a memory caching as
                // adjacent edges can be quickly computed
  }
  /*!
   * 			Incidence relations with faces
   * \param	edge	The involving edge
   * \return	An iterator range on the incident faces (AIFFace pointer) of the
   * involved edge.
   */
  static boost::iterator_range< incident_face_iterator >
  incident_faces(edge_descriptor edge)
  {
    if(edge == null_edge())
      throw std::invalid_argument("AIFTopologyHelpers::incident_faces -> no "
                                  "incident relation exists for a null edge.");

    return edge->GetIncidentFaces();
  }
  /*!
   * 			Function determining the position of a vertex in an edge
   * \param	edge	The involving edge
   * \param	vertex	The involving vertex
   * \return	The position of the vertex in the involving edge (FIRST or
   * SECOND).
   */
  static vertex_pos vertex_position(edge_descriptor edge,
                                    vertex_descriptor vertex)
  {
    if(edge->get_first_vertex() == vertex)
      return vertex_pos::FIRST;
    else if(edge->get_second_vertex() == vertex)
      return vertex_pos::SECOND;
    else
      throw std::invalid_argument(
          "AIFTopologyHelpers::vertex_position -> given vertex not found.");
  }
  /*!
   * 			Opposite vertex in an edge
   * \param	edge	The involving edge
   * \param	vertex	The involving vertex
   * \return	The first vertex of \c edge if \c vertex is the second, otherwise
   * the second vertex.
   */
  static vertex_descriptor opposite_vertex(edge_descriptor edge,
                                           vertex_descriptor vertex)
  {
    if(edge->get_first_vertex() == vertex)
      return edge->get_second_vertex();
    else if(edge->get_second_vertex() == vertex)
      return edge->get_first_vertex();
    else
      throw std::invalid_argument(
          "AIFTopologyHelpers::opposite_vertex -> given vertex not found.");
  }
  /*!
   * Vertices swap (first become second and second become first)
   *
   * \param  edge  The involving edge
   */
  static void swap_vertices(edge_descriptor edge)
  { // No std::swap cause get_first_vertex is a const function
    vertex_descriptor v_temp = edge->get_first_vertex();
    edge->set_first_vertex(edge->get_second_vertex());
    edge->set_second_vertex(v_temp);
  }
  /*!
   * Edge and face link.
   * An incidence relation is created between the edge and the face
   * by adding the edge in the incident edges in the edgecontainer of the face,
   * and adding the face in the incident faces container of the edge.
   *
   * \param  edge  The involving edge
   * \param  face   The involving face
   *
   * \warning  Be careful, this function can link any edge and face together
   *           but the involving edge's vertices have to be coherent,
   *           otherwise topological incoherences will be introduce.
   */
  static void link_edge_and_face(edge_descriptor edge, face_descriptor face)
  {
    bool face_has_that_incident_edge = are_incident(face, edge);
    bool edge_has_that_incident_face = are_incident(edge, face);

    if(!face_has_that_incident_edge || !edge_has_that_incident_face)
    {
      if(!face_has_that_incident_edge)
        add_edge_to_face(face, edge);
      if(!edge_has_that_incident_face)
        add_face_to_edge(edge, face);
    }
    else
      throw std::invalid_argument(
          "AIFTopologyHelpers::link_edge_and_face -> incident relation already "
          "exists between the two elements.");
  }

  /*!
   * Edge and face link.
   * An incidence relation is created between the edge and the face
   * by adding the edge in the incident edges in the edgecontainer of the face,
   * and adding the face in the incident faces container of the edge.
   *
   * \param  edge       The involving edge
   * \param  face       The involving face
   * \param  prev_edge  The existing edge around face after which the edge
   *                    "edge" must be inserted
   *
   * \warning  Be careful, this function can link any edge and face together
   *           but the involving edge's vertices have to be coherent,
   *           otherwise topological incoherences will be introduced.
   */
  static void
  link_edge_and_face_around_face_after_edge(edge_descriptor edge,
                                            face_descriptor face,
                                            edge_descriptor prev_edge)
  {
    bool face_has_that_incident_edge = are_incident(face, edge);
    bool edge_has_that_incident_face = are_incident(edge, face);

    if(!face_has_that_incident_edge || !edge_has_that_incident_face)
    {
      if(!face_has_that_incident_edge)
        add_edge_to_face_after_edge(face, prev_edge, edge);
      if(!edge_has_that_incident_face)
        add_face_to_edge(edge, face);
    }
    else
      throw std::invalid_argument(
          "AIFTopologyHelpers::link_edge_and_face_around_face_after_edge -> incident relation already "
          "exists between the two elements.");
  }

  /*!
   * Edges and face link.
   * An incidence relation is created between the edges and the face
   * by sequentially adding the edges in the incident edges in the
   * edgecontainer of the face, and adding the face in the incident faces
   * container of the edge.
   *
   * \param  edge       The involving edges
   * \param  face       The involving face
   * \param  prev_edge  The existing edge around face after which the
   *                    sequence of all edges "edges" must be inserted
   *
   * \warning  Be careful, this function can link any edge and face together
   *           but the involving edge's vertices have to be coherent,
   *           otherwise topological incoherences will be introduced.
   */
  static void link_edges_and_face_around_face_after_edge(
      std::vector< edge_descriptor > &edges,
      face_descriptor face,
      edge_descriptor prev_edge)
  {
    edge_descriptor current_edge;
    std::vector< edge_descriptor >::iterator it(edges.begin()),
        ite(edges.end());
    for(; it != ite; ++it)
    {
      edge_descriptor current_edge = *it;
      link_edge_and_face_around_face_after_edge(current_edge, face, prev_edge);
      prev_edge = current_edge;
    }
  }


  /*!
   * 			Edge and face unlink.
   * \param	edge		The involving edge
   * \param	face		The involving face
   */
  static void unlink_edge_and_face(edge_descriptor edge, face_descriptor face)
  {
    bool face_has_that_incident_edge = are_incident(face, edge);
    bool edge_has_that_incident_face = are_incident(edge, face);

    if(!face_has_that_incident_edge && !edge_has_that_incident_face)
      throw std::invalid_argument(
          "AIFTopologyHelpers::unlink_edge_and_face -> no incident relation "
          "existing between the two elements.");

    if(edge_has_that_incident_face)
    {
      remove_face_from_edge(edge, face);
    }

    if(face_has_that_incident_edge)
    {
      remove_edge_from_face(face, edge);
    }
  }
  /*!
   * 			Edges and face unlink.
   * \param	edges		The involving edges
   * \param	face		The involving face
   */
  static void unlink_edges_and_face(std::vector< edge_descriptor > edges,
                                    face_descriptor face)
  {
    std::vector< edge_descriptor >::iterator it(edges.begin()),
        ite(edges.end());
    for(; it != ite; ++it)
      unlink_edge_and_face(*it, face);
  }

  /*!
   * Function determining if the argument edge is a 2-manifold edge
   * considering only the topology (not the geometry)
   *
   * \param  edge  The involving edge
   *
   * \return  true if the argument edge is a locally 2-manifold edge, false
   *          otherwise.
   */
  static bool is_2_manifold_edge(edge_descriptor edge)
  {
    if(is_degenerated_edge(edge))
      return false;

    return is_2_manifold_vertex(edge->get_first_vertex()) &&
           is_2_manifold_vertex(edge->get_second_vertex());
  }

  /*!
  * 			Sort mesh edges
  * \param	mesh	The involving mesh
  * \param	cmp   	The involving comparator object used for the sorting
  * \return	An iterator range on the edges (AIFEdge pointer) of the involve mesh.
  */
  template<typename ComparatorType>
  static boost::iterator_range<edge_container::const_iterator> sort_edges(ptr_mesh mesh, ComparatorType cmp)
  {
    std::sort(mesh->GetEdges().begin(), mesh->GetEdges().end(), cmp);
    return mesh->GetEdges();
  }

  /*!
  * 			Sort mesh edges
  * \param	mesh	The involving mesh
  * \param	cmp   	The involving comparator object used for the sorting
  * \return	An iterator range on the edges (AIFEdge pointer) of the involve mesh.
  */
  template<typename ComparatorType>
  static boost::iterator_range<edge_container::const_iterator> sort_edges(smart_ptr_mesh mesh, ComparatorType cmp)
  {
    return sort_edges<ComparatorType>(mesh.get(), cmp);
  }

  /*!
  * 			Sort mesh edges
  * \param	mesh	The involving mesh
  * \param	cmp   	The involving comparator object used for the sorting
  * \return	An iterator range on the edges (AIFEdge pointer) of the involve mesh.
  */
  template<typename ComparatorType>
  static boost::iterator_range<edge_container::const_iterator> sort_edges(ref_mesh mesh, ComparatorType cmp)
  {
    return sort_edges<ComparatorType>(&mesh, cmp);
  }
  //////////////////////////////////	Face AIFTopologyHelpers
  /////////////////////////////////////
  /*!
   * 			Face degree
   * \param	face	The involving face
   * \return	The number of incident edges of the argument face.
   */
  static size_type degree(face_descriptor face) { return face->GetDegree(); }

  /*!
   * Function determining if the argument face is isolated (without any
   * adjacent face + its vertices are degree 2 vertices)
   *
   * \param  face  The involving face
   *
   * \return  true if the argument face is isolated, false otherwise.
   */
  static bool is_isolated_face(face_descriptor face)
  {
    auto edges_range = incident_edges(face);
    auto it = edges_range.begin();
    for(; it != edges_range.end(); ++it)
    {
      if((degree((*it)->get_first_vertex()) > 2) ||
         (degree((*it)->get_second_vertex()) > 2))
        return false;
    }
    return true;
  }
  /*!
   * 			Function determining if the argument face is degenerated
   * \param	face	The involving face
   * \return	true if the argument face is null, or its degree is strictly less
   * than 3, or if it has a degenerated edge among its incident edges, false
   * otherwise.
   */
  static bool is_degenerated_face(face_descriptor face)
  {
    if((face == null_face()) || // a null face is degenerated
       (degree(face) <
        3) // a face with degree strictly less than 3 is degenerated
    )
      return true;
    auto eRange = incident_edges(face);
    for(auto eIt = eRange.begin(); eIt != eRange.end(); ++eIt)
    {
      if(is_degenerated_edge(*eIt))
        return true;
    }
    return false;
  }

  /*!
   * Function determining if the argument face is a 2-manifold face
   * considering only the topology (not the geometry)
   *
   * \param  face  The involving face
   *
   * \return  true if the argument face is a locally 2-manifold face, false
   *          otherwise.
   */
  static bool is_2_manifold_face(face_descriptor face)
  {
    if(is_degenerated_face(face))
      return false;

    auto edges_range = incident_edges(face);
    auto it = edges_range.begin();
    for(; it != edges_range.end(); ++it)
    {
      if(!is_2_manifold_edge(*it))
        return false;
    }
    return true;
  }

  /*!
  * 			Function determining if the argument face is dangling (not isolated and without any adjacent face sharing and edge)
  * \param	face	The involving face
  * \return	true if the argument face is dangling, false otherwise.
  */
  static bool is_dangling_face(face_descriptor face)
  {
    bool at_least_one_vertex_with_degree_greater_than_2 = false;
    auto edges_range = incident_edges(face);
    auto it = edges_range.begin();
    for (; it != edges_range.end(); ++it) {
      if (degree(*it) > 1)
        return false;
      if (!at_least_one_vertex_with_degree_greater_than_2 && 
          ((degree((*it)->get_first_vertex()) > 2) || 
           (degree((*it)->get_second_vertex()) > 2)))
        at_least_one_vertex_with_degree_greater_than_2 = true;
    }
    return at_least_one_vertex_with_degree_greater_than_2;
  }

  /*!
  * 			Function determining if the argument face has only incident edges with degree equals to 1
  * \param	face	The involving face
  * \return	true if the argument face has only incident edges with degree equals to 1.
  */
  static bool face_has_only_incident_edges_with_degree_1(face_descriptor face)
  {
    auto edges_range = incident_edges(face);
    auto it = edges_range.begin();
    for (; it != edges_range.end(); ++it) {
      if (degree(*it) > 1)
        return false;
    }
    return true;
  }

  /*!
  * 			Function determining if the argument face has adjacent faces only along one incident edge
  * \param	face	The involving face
  * \param	edge	The involving edge
  * \return	true if the argument face has only the edge edge with degree greater than 1.
  */
  static bool face_has_only_that_edge_with_degree_greater_than_1(face_descriptor face, edge_descriptor edge)
  {
    auto edges_range = incident_edges(face);
    auto it = edges_range.begin();
    if (!are_incident(face, edge))
      return false;
    if (degree(edge) <= 1)
      return false;
    for (; it != edges_range.end(); ++it) {
      if ((degree(*it) > 1) && (*it != edge))
        return false;
    }
    return true;
  }

  /*!
   * \brief  Function determining if the argument face shares an incidence
   *         relation with the argument vertex.
   *
   * \param  face    The involving face
   * \param  vertex  The involving vertex
   *
   * \return true if the arguments share an incidence relation, false
   *         otherwise.
   */
  static bool are_incident(face_descriptor face, vertex_descriptor vertex)
  {
	if(face == null_face())
      return false;		  
    if(face->m_Is_Incident_PtrVertices_Computed) // optimization by using
                                                 // topological caching
    {
      return (std::find(face->m_Incident_PtrVertices.begin(),
                        face->m_Incident_PtrVertices.end(),
                        vertex) != face->m_Incident_PtrVertices.end());
    }
    else
    {
      auto edges_range = incident_edges(face);
      auto it = edges_range.begin();
      for(; it != edges_range.end(); ++it)
      {
        if(are_incident(*it, vertex))
          return true;
      }
    }
    return false;
  }

  /*!
   * \brief  Function determining if the argument vertex shares an incidence
   *         relation with the argument face.
   *
   * \param  vertex  The involving vertex
   * \param  face    The involving face
   *
   * \return true if the arguments share an incidence relation, false
   *         otherwise.
   */
  static bool are_incident(vertex_descriptor vertex, face_descriptor face)
  {
	if(vertex == null_vertex())
      return false;		  	  
    if(vertex->m_Incident_PtrFaces_Computed) // optimization by using
                                             // topological caching
    {
      return (std::find(vertex->m_Incident_PtrFaces.begin(),
                        vertex->m_Incident_PtrFaces.end(),
                        face) != vertex->m_Incident_PtrFaces.end());
    }
    else
    {
      auto edges_range = incident_edges(vertex);
      auto it = edges_range.begin();
      for(; it != edges_range.end(); ++it)
      {
        if(are_incident(*it, face))
          return true;
      }
    }
    return false;
  }
  /*!
   * \brief  Function determining if the argument face shares an incidence
   *         relation with the argument edge.
   *
   * \param  face  The involving face
   * \param  edge  The involving edge
   *
   * \return  true if the arguments share an incidence relation, false
   *          otherwise.
   *
   * \warning are_incident(face, edge) does not give necessary the
   *          same result as are_incident(edge, face) when topological
   *          modifications are in progress.
   */
  static bool are_incident(face_descriptor face, edge_descriptor edge)
  {
	if(face == null_face())
      return false;		  
    auto edges_range = incident_edges(face);
    auto it = edges_range.begin();
    for(; it != edges_range.end(); ++it)
    {
      if(*it == edge)
        return true;
    }
    return false;
  }
  /*!
   * Function determining if the arguments share an adjacent relation
   * (if they share an edge)
   *
   * \param  face1  The first involving face
   * \param  face2  The first involving face
   *
   * \return  true if the arguments share an adjacent relation, false
   *          otherwise (a face is not adjacent to itself).
   */
  static bool are_adjacent(face_descriptor face1, face_descriptor face2)
  {
    if(face1 == face2)
      return false;

    boost::iterator_range< incident_edge_iterator > edges_range =
        incident_edges(face1);

    for(incident_edge_iterator it = edges_range.begin();
        it != edges_range.end();
        ++it)
    {
      if(are_incident(*it, face2))
        return true;
    }
    return false;
  }
  /*!
   * 			Incidence relations with vertices
   * \param	face	The involving face
   * \return	An iterator range on the incident vertices (AIFVertex pointer) of
   * the involved face.
   */
  static boost::iterator_range< vertex_container_in_face::const_iterator >
  incident_vertices(face_descriptor face)
  {
    if(face == null_face())
      throw std::invalid_argument("AIFTopologyHelpers::incident_vertices -> no "
                                  "incident relation exists for a null face.");

    return face->GetIncidentVertices();
  }
  /*!
   * 			Incidence relations with edges
   * \param	face	The involving face
   * \return	An iterator range on the incident edges (AIFEdge pointer) of the
   * involved face.
   */
  static boost::iterator_range< incident_edge_iterator >
  incident_edges(face_descriptor face)
  {
    if(face == null_face())
      throw std::invalid_argument("AIFTopologyHelpers::incident_edges -> no "
                                  "incident relation exists for a null face.");

    return face->GetIncidentEdges();
  }
  /*!
   * 			Adjacent relations with faces
   * \param	face	The involving face
   * \return	An iterator range on the adjacent faces (AIFFace pointer) of the
   * involved face. When one_neighborhood_sharing_edge is true, the adjacent
   * faces are those sharing an edge with involved face. Else the adjacent faces
   * are those sharing a vertex with involved face.
   */
  static std::vector< face_descriptor >
  adjacent_faces(face_descriptor face,
                 bool one_neighborhood_sharing_edge = true)
  {
    std::set< face_descriptor > adFaces;
    if(one_neighborhood_sharing_edge)
    {
      auto eRange = incident_edges(face);
      for(auto eIt = eRange.begin(); eIt != eRange.end(); ++eIt)
      {
        auto fRange = incident_faces(*eIt);
        std::vector< face_descriptor > incidentFaces(fRange.begin(),
                                                     fRange.end());
        incidentFaces.erase(
            std::find(incidentFaces.begin(), incidentFaces.end(), face));
        adFaces.insert(incidentFaces.begin(), incidentFaces.end());
      }
    }
    else
    {
      auto vRange = incident_vertices(face);
      for(auto vIt = vRange.begin(); vIt != vRange.end(); ++vIt)
      {
        auto fRange = incident_faces(*vIt);
        std::vector< face_descriptor > incidentFaces(fRange.begin(),
                                                     fRange.end());
        incidentFaces.erase(
            std::find(incidentFaces.begin(), incidentFaces.end(), face));
        adFaces.insert(incidentFaces.begin(), incidentFaces.end());
      }
    }
    return std::vector< face_descriptor >(
        adFaces.begin(),
        adFaces.end()); // note that it is a non-sens to manage a memory caching
                        // as adjacent faces can be quickly computed
  }

  /*!
  * 			Reverse the incident edge order of the given face.
  * \param	face	The involved face
  */
  static void reverse_face_orientation(face_descriptor face)
  {
    auto edgeRange = incident_edges(face);
    std::reverse(edgeRange.begin(), edgeRange.end());
    face->clear_vertex_incidency();
  }

  /*!
  * 			First common face between two edges
  * \param	edge1	The first involving edge
  * \param	edge2	The second involving edge
  * \return	The first common face between edge1 and edge2.
  */
  static face_descriptor common_face(edge_descriptor edge1, edge_descriptor edge2)
  {
    typedef face_container_in_edge::const_iterator it_type;

    boost::iterator_range<it_type> faces_range = incident_faces(edge1);

    for (it_type it = faces_range.begin(); it != faces_range.end(); ++it) {
      if (are_incident(*it, edge2))
        return *it;
    }
    return null_face();
  }
  //////////////////////////////////	Mesh AIFTopologyHelpers
  /////////////////////////////////////
  /*!
   * 			Vertices count getter
   * \param	mesh	The involving mesh
   * \return	The number of vertices of the involved mesh.
   */
  static size_type num_vertices(ptr_cmesh mesh)
  {
    return mesh->GetNumberOfVertices();
  }
  /*!
   * 			Vertices count getter
   * \param	mesh	The involving mesh
   * \return	The number of vertices of the involved mesh.
   */
  static size_type num_vertices(smart_ptr_cmesh mesh)
  {
    return num_vertices(mesh.get());
  }
  /*!
   * 			Vertices count getter
   * \param	mesh	The involving mesh
   * \return	The number of vertices of the involved mesh.
   * \warning  const qualifier on reference has no effect for MAC OS
   */
  static size_type num_vertices(ref_cmesh mesh) { return num_vertices(&mesh); }

  /*!
   * 			Edges count getter
   * \param	mesh	The involving mesh
   * \return	The number of edges of the involved mesh.
   */
  static size_type num_edges(ptr_cmesh mesh)
  {
    return mesh->GetNumberOfEdges();
  }
  /*!
   * 			Edges count getter
   * \param	mesh	The involving mesh
   * \return	The number of edges of the involved mesh.
   */
  static size_type num_edges(smart_ptr_cmesh mesh)
  {
    return num_edges(mesh.get());
  }
  /*!
   * 			Edges count getter
   * \param	mesh	The involving mesh
   * \return	The number of edges of the involved mesh.
   */
  static size_type num_edges(ref_cmesh mesh) { return num_edges(&mesh); }


  /*!
   * 			Faces count getter
   * \param	mesh	The involving mesh
   * \return	The number of faces of the involved mesh.
   */
  static size_type num_faces(ptr_cmesh mesh)
  {
    return mesh->GetNumberOfFaces();
  }
  /*!
   * 			Faces count getter
   * \param	mesh	The involving mesh
   * \return	The number of faces of the involved mesh.
   */
  static size_type num_faces(smart_ptr_cmesh mesh)
  {
    return num_faces(mesh.get());
  }
  /*!
   * 			Faces count getter
   * \param	mesh	The involving mesh
   * \return	The number of faces of the involved mesh.
   */
  static size_type num_faces(ref_cmesh mesh) { return num_faces(&mesh); }

  /*!
   * 			Mesh vertices getter
   * \param	mesh	The involving mesh
   * \return	An iterator range on the vertices (AIFVertex pointer) of the
   * involved mesh.
   */
  static boost::iterator_range< vertex_container::const_iterator >
  vertices(ptr_cmesh mesh)
  {
    return mesh->GetVertices();
  }
  /*!
   * 			Mesh vertices getter
   * \param	mesh	The involving mesh
   * \return	An iterator range on the vertices (AIFVertex pointer) of the
   * involved mesh.
   */
  static boost::iterator_range< vertex_container::const_iterator >
  vertices(smart_ptr_cmesh mesh)
  {
    return vertices(mesh.get());
  }
  /*!
   * 			Mesh vertices getter
   * \param	mesh	The involving mesh
   * \return	An iterator range on the vertices (AIFVertex pointer) of the
   * involved mesh.
   */
  static boost::iterator_range< vertex_container::const_iterator >
  vertices(ref_cmesh mesh)
  {
    return vertices(&mesh);
  }

  /*!
   * 			Mesh edges getter
   * \param	mesh	The involving mesh
   * \return	An iterator range on the edges (AIFEdge pointer) of the involved
   * mesh.
   */
  static boost::iterator_range< edge_container::const_iterator >
  edges(ptr_cmesh mesh)
  {
    return mesh->GetEdges();
  }
  /*!
   * 			Mesh edges getter
   * \param	mesh	The involving mesh
   * \return	An iterator range on the edges (AIFEdge pointer) of the involved
   * mesh.
   */
  static boost::iterator_range< edge_container::const_iterator >
  edges(smart_ptr_cmesh mesh)
  {
    return edges(mesh.get());
  }
  /*!
   * 			Mesh edges getter
   * \param	mesh	The involving mesh
   * \return	An iterator range on the edges (AIFEdge pointer) of the involved
   * mesh.
   */
  static boost::iterator_range< edge_container::const_iterator >
  edges(ref_cmesh mesh)
  {
    return edges(&mesh);
  }

  /*!
   * Function determining if the argument mesh is a 2-manifold
   * considering only the topology (not the geometry).
   *
   * \param  vertex  The involving mesh
   *
   * \return  true if the argument mesh is a 2-manifold mesh,
   *          false otherwise.
   */
  static bool is_2_manifold_mesh(ptr_cmesh mesh)
  {
    auto verticesRange = mesh->GetVertices();
    auto vertexIter = verticesRange.begin(), vertexIterE = verticesRange.end();
    for(; vertexIter != vertexIterE; ++vertexIter)
    {
      if(!is_2_manifold_vertex(*vertexIter))
        return false;
    }
    return true;
  }

  /*!
   * Function determining if the argument mesh is a 2-manifold
   * considering only the topology (not the geometry).
   *
   * \param  vertex  The involving mesh
   *
   * \return  true if the argument mesh is a 2-manifold mesh, false otherwise.
   */
  static bool is_2_manifold_mesh(smart_ptr_cmesh mesh)
  {
    return edges(mesh.get());
  }

  /*!
   * Function determining if the argument mesh is a 2-manifold
   * considering only the topology (not the geometry).
   *
   * \param  vertex  The involving mesh
   *
   * \return  true if the argument mesh is a 2-manifold mesh,
   *          false otherwise.
   */
  static bool is_2_manifold_mesh(ref_cmesh mesh)
  {
    return is_2_manifold_mesh(&mesh);
  }


  /*!
   *\return	 true if the border edges are in normalized representation, which
   *is when enumerating all edges with the iterator : The non - border edges
   *precede the border edges.
   */
  static bool normalized_border_is_valid(ptr_cmesh mesh)
  {
    auto edgesRange = mesh->GetEdges();
    auto edgeIter = edgesRange.begin(), edgeIterE = edgesRange.end();
    bool meetFirstBorderEdge = false;
    for(; edgeIter != edgeIterE; ++edgeIter)
    {
      if(is_surface_border_edge(*edgeIter))
      {
        if(!meetFirstBorderEdge)
          meetFirstBorderEdge = true;
      }
      else
      {
        if(meetFirstBorderEdge)
          return false;
      }
    }
    return true;
  }
  /*!
   *\return	 true if the border edges are in normalized representation, which
   *is when enumerating all edges with the iterator : The non - border edges
   *precede the border edges.
   */
  static bool normalized_border_is_valid(smart_ptr_cmesh mesh)
  {
    return normalized_border_is_valid(mesh.get());
  }
  /*!
   *\return	 true if the border edges are in normalized representation, which
   *is when enumerating all edges with the iterator : The non - border edges
   *precede the border edges.
   */
  static bool normalized_border_is_valid(ref_cmesh mesh)
  {
    return normalized_border_is_valid(&mesh);
  }

  /*!
   * \return the number of border edges.
   */
  static size_type size_of_border_edges(ptr_cmesh mesh)
  {
    auto edgesRange = mesh->GetEdges();
    auto edgeIter = edgesRange.begin(), edgeIterE = edgesRange.end();
    size_type cpt = 0;
    for(; edgeIter != edgeIterE; ++edgeIter)
    {
      if(is_surface_border_edge(*edgeIter))
      {
        ++cpt;
      }
    }
    return cpt;
  }
  /*!
   * \return the number of border edges.
   */
  static size_type size_of_border_edges(smart_ptr_cmesh mesh)
  {
    return size_of_border_edges(mesh.get());
  }
  /*!
   * \return the number of border edges.
   */
  static size_type size_of_border_edges(ref_cmesh mesh)
  {
    return size_of_border_edges(&mesh);
  }

  /*!
   * Sorts edges such that the non-border edges precede the border
   * edges.
   *
   * \param  mesh  The involving mesh
   */
  static void normalize_border(ptr_mesh mesh)
  {
    auto edgesRange = mesh->GetEdges();
    auto edgeIter = edgesRange.begin(), edgeIterE = edgesRange.end();
    --edgeIterE;
    while(edgeIter != edgeIterE)
    {
      if(is_surface_border_edge(*edgeIter))
      {
        if(!is_surface_border_edge(*edgeIterE))
        {
          swap(*edgeIter, *edgeIterE);
          ++edgeIter;
        }
        --edgeIterE;
      }
      else
        ++edgeIter;
    }
    assert(normalized_border_is_valid(mesh));
  }
  /*!
   * Sorts edges such that the non-border edges precede the border
   * edges.
   *
   * \param  mesh  The involving mesh
   */
  static void normalize_border(smart_ptr_mesh mesh)
  {
    normalize_border(mesh.get());
  }
  /*!
   * Sorts edges such that the non-border edges precede the border
   * edges.
   *
   * \param  mesh  The involving mesh
   */
  static void normalize_border(ref_mesh mesh) { normalize_border(&mesh); }

  /*!
   * 			Mesh border edges getter
   * \param	mesh	The involving mesh
   * \return	An iterator range on the border edges (AIFEdge pointer) of the
   * involved mesh.
   */
  static boost::iterator_range< edge_container::const_iterator >
  border_edges(ptr_cmesh mesh)
  {
    if(!normalized_border_is_valid(mesh))
      throw std::invalid_argument(
          "AIFTopologyHelpers::border_edges(m) -> mesh' edges have not been "
          "normalized, thus cannot deduced an iterator range over border "
          "edges.");

    auto edgesRange = mesh->GetEdges();
    auto edgeIter = edgesRange.begin(), edgeIterE = edgesRange.end();
    while((edgeIter != edgeIterE) && (!is_surface_border_edge(*edgeIter)))
      ++edgeIter;
    return boost::make_iterator_range(edgeIter, edgeIterE);
  }
  /*!
   * 			Mesh border edges getter
   * \param	mesh	The involving mesh
   * \return	An iterator range on the border edges (AIFEdge pointer) of the
   * involved mesh.
   */
  static boost::iterator_range< edge_container::const_iterator >
  border_edges(smart_ptr_cmesh mesh)
  {
    return border_edges(mesh.get());
  }
  /*!
   * 			Mesh border edges getter
   * \param	mesh	The involving mesh
   * \return	An iterator range on the border edges (AIFEdge pointer) of the
   * involved mesh.
   */
  static boost::iterator_range< edge_container::const_iterator >
  border_edges(ref_cmesh mesh)
  {
    return border_edges(&mesh);
  }

  /*!
   * 			Mesh faces getter
   * \param	mesh	The involving mesh
   * \return	An iterator range on the faces (AIFFace pointer) of the involved
   * mesh.
   */
  static boost::iterator_range< face_container::const_iterator >
  faces(ptr_cmesh mesh)
  {
    return mesh->GetFaces();
  }
  /*!
   * 			Mesh faces getter
   * \param	mesh	The involving mesh
   * \return	An iterator range on the faces (AIFFace pointer) of the involved
   * mesh.
   */
  static boost::iterator_range< face_container::const_iterator >
  faces(smart_ptr_cmesh mesh)
  {
    return faces(mesh.get());
  }
  /*!
   * 			Mesh faces getter
   * \param	mesh	The involving mesh
   * \return	An iterator range on the faces (AIFFace pointer) of the involved
   * mesh.
   */
  static boost::iterator_range< face_container::const_iterator >
  faces(ref_cmesh mesh)
  {
    return faces(&mesh);
  }

  /*!
  * 			Generic function to know if an iterator range contains a
  *       degenerated face.
  * \param	first	The iterator pointing on the range beginning.
  * \param	last	The iterator pointing on the range end.
  * \return	True if the iterator range contains a degenerated face, false
            otherwise.
  */
  template< typename InputIt >
  static bool contains_a_degenerated_face(InputIt first, InputIt last)
  {
    for(; first != last; ++first)
      if(is_degenerated_face(*first))
        return true;

    return false;
  }

  /*!
   * 			Mesh vertex inserter
   * \param	mesh	The involving mesh
   * \return	The vertex newly added to the involved mesh
   */
  static vertex_descriptor add_vertex(ptr_mesh mesh)
  {
    vertex_descriptor vertex = vertex_type::New();
    mesh->InsertIsolatedVertex(vertex);
    return vertex;
  }
  /*!
   * 			Mesh vertex inserter
   * \param	mesh	The involving mesh
   * \return	The vertex newly added to the involved mesh
   */
  static vertex_descriptor add_vertex(smart_ptr_mesh mesh)
  {
    return add_vertex(mesh.get());
  }
  /*!
   * 			Mesh vertex inserter
   * \param	mesh	The involving mesh
   * \return	The vertex newly added to the involved mesh
   */
  static vertex_descriptor add_vertex(ref_mesh mesh)
  {
    return add_vertex(&mesh);
  }

  /*!
   * 			Mesh vertex inserter
   * \param	vertex	The isolated vertex to insert
   * \param	mesh	The involving mesh
   * \return	void
   */
  static void add_vertex(vertex_descriptor vertex, ptr_mesh mesh)
  {
    mesh->InsertIsolatedVertex(vertex);
  }

  /*!
   * 			Mesh vertex inserter
   * \param	vertex	The isolated vertex to insert
   * \param	mesh	The involving mesh
   * \return	void
   */
  static void add_vertex(vertex_descriptor vertex, smart_ptr_mesh mesh)
  {
    add_vertex(vertex, mesh.get());
  }

  /*!
   * 			Mesh vertex inserter
   * \param	vertex	The isolated vertex to insert
   * \param	mesh	The involving mesh
   * \return	void
   */
  static void add_vertex(vertex_descriptor vertex, ref_mesh mesh)
  {
    add_vertex(vertex, &mesh);
  }

  /*!
   * 			Mesh edge inserter
   * \param	mesh	The involving mesh
   * \return	The edge newly added to the involved mesh
   */
  static edge_descriptor add_edge(ptr_mesh mesh)
  {
    edge_descriptor edge = edge_type::New();
    mesh->InsertIsolatedEdge(edge);
    return edge;
  }
  /*!
   * 			Mesh edge inserter
   * \param	mesh	The involving mesh
   * \return	The edge newly added to the involved mesh
   */
  static edge_descriptor add_edge(smart_ptr_mesh mesh)
  {
    return add_edge(mesh.get());
  }
  /*!
   * 			Mesh edge inserter
   * \param	mesh	The involving mesh
   * \return	The edge newly added to the involved mesh
   */
  static edge_descriptor add_edge(ref_mesh mesh) { return add_edge(&mesh); }

  /*!
   * 			Mesh edge inserter
   * \param	edge	The isolated edge to insert
   * \param	mesh	The involving mesh
   * \return	void
   */
  static void add_edge(edge_descriptor edge, ptr_mesh mesh)
  {
    mesh->InsertIsolatedEdge(edge);
  }

  /*!
   * 			Mesh edge inserter
   * \param	edge	The isolated edge to insert
   * \param	mesh	The involving mesh
   * \return	void
   */
  static void add_edge(edge_descriptor edge, smart_ptr_mesh mesh)
  {
    add_edge(edge, mesh.get());
  }

  /*!
   * 			Mesh edge inserter
   * \param	edge	The isolated edge to insert
   * \param	mesh	The involving mesh
   * \return	void
   */
  static void add_edge(edge_descriptor edge, ref_mesh mesh)
  {
    add_edge(edge, &mesh);
  }

  /*!
   * 			Mesh edge inserter
   * \param	v	A vertex incident to the edge to insert
   * \param	u	Another vertex incident to the edge to insert
   * \param	mesh	The involving mesh
   * \return	std::pair<edge_descriptor,bool>; if the edge (u,v) is already in
   * the graph, then the bool flag returned is false and the returned edge
   * descriptor points to the already existing edge
   */
  static std::pair< edge_descriptor, bool >
  add_edge(vertex_descriptor v, vertex_descriptor u, ptr_mesh mesh)
  {
    bool is_edge_insertion_done = false;
    edge_descriptor edge = common_edge(v, u);
    if(edge == null_edge())
    {
      edge = add_edge(mesh); // does the mesh->InsertIsolatedEdge(edge);
      link_vertex_and_edge(v, edge, vertex_pos::FIRST);
      link_vertex_and_edge(u, edge, vertex_pos::SECOND);
      is_edge_insertion_done = true;
    }

    return std::pair< edge_descriptor, bool >(edge, is_edge_insertion_done);
  }

  /*!
   * 			Mesh edge inserter
   * \param	v	A vertex incident to the edge to insert
   * \param	u	Another vertex incident to the edge to insert
   * \param	mesh	The involving mesh
   * \return	std::pair<edge_descriptor,bool>; if the edge (u,v) is already in
   * the graph, then the bool flag returned is false and the returned edge
   * descriptor points to the already existing edge
   */
  static std::pair< edge_descriptor, bool >
  add_edge(vertex_descriptor v, vertex_descriptor u, smart_ptr_mesh mesh)
  {
    return add_edge(v, u, mesh.get());
  }

  /*!
   * 			Mesh edge inserter
   * \param	v	A vertex incident to the edge to insert
   * \param	u	Another vertex incident to the edge to insert
   * \param	mesh	The involving mesh
   * \return	std::pair<edge_descriptor,bool>; if the edge (u,v) is already in
   * the graph, then the bool flag returned is false and the returned edge
   * descriptor points to the already existing edge
   */
  static std::pair< edge_descriptor, bool >
  add_edge(vertex_descriptor v, vertex_descriptor u, ref_mesh mesh)
  {
    return add_edge(v, u, &mesh);
  }

  /*!
   * 			Mesh face inserter
   * \param	mesh	The involving mesh
   * \return	The face newly added to the involved mesh
   */
  static face_descriptor add_face(ptr_mesh mesh)
  {
    face_descriptor face = face_type::New();
    mesh->InsertIsolatedFace(face);
    return face;
  }
  /*!
   * 			Mesh face inserter
   * \param	mesh	The involving mesh
   * \return	The face newly added to the involved mesh
   */
  static face_descriptor add_face(smart_ptr_mesh mesh)
  {
    return add_face(mesh.get());
  }
  /*!
   * 			Mesh face inserter
   * \param	mesh	The involving mesh
   * \return	The face newly added to the involved mesh
   */
  static face_descriptor add_face(ref_mesh mesh) { return add_face(&mesh); }
  /*!
   * 			Mesh face inserter
   * \param	face	The isolated face to insert
   * \param	mesh	The involving mesh
   * \return	The face newly added to the involved mesh
   */
  static void add_face(face_descriptor face, ptr_mesh mesh)
  {
    mesh->InsertIsolatedFace(face);
  }
  /*!
   * 			Mesh face inserter
   * \param	face	The isolated face to insert
   * \param	mesh	The involving mesh
   * \return	The face newly added to the involved mesh
   */
  static void add_face(face_descriptor face, smart_ptr_mesh mesh)
  {
    add_face(face, mesh.get());
  }
  /*!
   * 			Mesh face inserter
   * \param	face	The isolated face to insert
   * \param	mesh	The involving mesh
   * \return	The face newly added to the involved mesh
   */
  static void add_face(face_descriptor face, ref_mesh mesh)
  {
    add_face(face, &mesh);
  }
  /////////////////////////////////
  /*!
   * Function to clear vertex'edge incidency. Equivalent to
   * the unlink_all_edges(vertex) function. The only interest of this function
   * is to have a function following CGAL naming conventions.
   *
   * \param  vertex  The involving vertex
   */
  static void clear_vertex(vertex_descriptor vertex)
  {
    if(vertex == null_vertex())
      throw std::invalid_argument(
          "AIFTopologyHelpers::clear_vertex(v) -> vertex is null, cannot clear "
          "its connectivity from mesh.");
    unlink_all_edges(vertex);
  }
  /*!
   * 			Function to clear vertex'edge incidency if present
   *       and then to remove the involved vertex from the mesh.
   * \param	vertex	The involving vertex
   * \param	mesh	The involving mesh
   */
  static void remove_vertex(vertex_descriptor vertex, ptr_mesh mesh)
  {
    if(vertex == null_vertex())
      throw std::invalid_argument(
          "AIFTopologyHelpers::remove_vertex(v, m) -> vertex is null, cannot "
          "suppress it from mesh.");
    clear_vertex(vertex); // in common mesh datastructures this check is not
                          // done beacause a call to clear_vertex is supposed to
                          // done just before a call to remove_vertex
    mesh->EraseIsolatedVertex(vertex);
  }
  /*!
   * 			Function to clear vertex'edge incidency if present
   *       and then to remove the involved vertex from the mesh.
   * \param	vertex	The involving vertex
   * \param	mesh	The involving mesh
   */
  static void remove_vertex(vertex_descriptor vertex, smart_ptr_mesh mesh)
  {
    remove_vertex(vertex, mesh.get());
  }
  /*!
   * 			Function to clear vertex'edge incidency if present
   *       and then to remove the involved vertex from the mesh.
   * \param	vertex	The involving vertex
   * \param	mesh	The involving mesh
   */
  static void remove_vertex(vertex_descriptor vertex, ref_mesh mesh)
  {
    remove_vertex(vertex, &mesh);
  }
  /*!
   * Generic function to remove vertices given by an iterator
   * range.
   *
   * \param  first  The iterator pointing on the range beginning.
   * \param  last   The iterator pointing on the range end.
   * \param  mesh   The involving mesh.
   */
  template< typename InputIt, typename MeshType >
  static void remove_vertices(InputIt first, InputIt last, MeshType& mesh)
  {
    for(; first != last; ++first)
      remove_vertex(*first, mesh);
  }
  /*!
   * 			Function to clear edge'face incidency and vertex incidency if
   * present and then to remove the involved edge from the mesh. Resulting
   *       isolated vertices are removed.
   * \param	edge	The involving edge.
   * \param	mesh	The involving mesh.
   */
  static void remove_edge(edge_descriptor edge, ptr_mesh mesh)
  {
    if(edge == null_edge())
      throw std::invalid_argument(
          "AIFTopologyHelpers::remove_edge(e, m) -> edge is null, cannot "
          "suppress it from mesh.");
    unlink_all_faces(edge);
    std::set< vertex_descriptor > verticesToRemove;
    if(edge->get_first_vertex() != null_vertex() &&
       edge->get_first_vertex()->GetDegree() == 1)
      verticesToRemove.insert(edge->get_first_vertex());
    if(edge->get_second_vertex() != null_vertex() &&
       edge->get_second_vertex()->GetDegree() == 1)
      verticesToRemove.insert(edge->get_second_vertex());
    unlink_all_vertices(edge);
    mesh->EraseIsolatedEdge(edge);

    std::set< vertex_descriptor >::iterator itv = verticesToRemove.begin(),
                                            itve = verticesToRemove.end();
    for(; itv != itve; ++itv)
    {
      mesh->EraseIsolatedVertex(*itv);
    }
  }
  /*!
   * 			Function to clear edge'face incidency and vertex incidency if
   * present and then to remove the involved edge from the mesh. Resulting
   *       isolated vertices are removed.
   * \param	edge	The involving edge.
   * \param	mesh	The involving mesh.
   */
  static void remove_edge(edge_descriptor edge, smart_ptr_mesh mesh)
  {
    remove_edge(edge, mesh.get());
  }
  /*!
   * 			Function to clear edge'face incidency and vertex incidency if
   * present and then to remove the involved edge from the mesh. Resulting
   *       isolated vertices are removed.
   * \param	edge	The involving edge.
   * \param	mesh	The involving mesh.
   */
  static void remove_edge(edge_descriptor edge, ref_mesh mesh)
  {
    remove_edge(edge, &mesh);
  }
  /*!
   * Generic function to remove edges given by an iterator
   * range.
   *
   * \param  first  The iterator pointing on the range beginning.
   * \param  last   The iterator pointing on the range end.
   * \param  mesh   The involving mesh.
   */
  template< typename InputIt, typename MeshType >
  static void remove_edges(InputIt first, InputIt last, MeshType& mesh)
  {
    for(; first != last; ++first)
      remove_edge(*first, mesh);
  }
  /*!
   * 			Function to clear edge'face incidency and vertex incidency if
   * present and then to remove the involved edge from the mesh. Resulting
   *       isolated vertices are removed.
   * \param	v	The first vertex of the involved edge.
   * \param	u	The second vertex of the involved edge.
   * \param	mesh	The involving mesh.
   */
  static void
  remove_edge(vertex_descriptor v, vertex_descriptor u, ptr_mesh mesh)
  {
    edge_descriptor edge = common_edge(v, u);
    remove_edge(edge, mesh); // to ensure same behavior with
                             // remove_edge(edge_descriptor edge, ptr_mesh mesh)
  }
  /*!
   * 			Function to clear edge'face incidency and vertex incidency if
   * present and then to remove the involved edge from the mesh. Resulting
   *       isolated vertices are removed.
   * \param	v	The first vertex of the involved edge.
   * \param	u	The second vertex of the involved edge.
   * \param	mesh	The involving mesh.
   */
  static void
  remove_edge(vertex_descriptor v, vertex_descriptor u, smart_ptr_mesh mesh)
  {
    remove_edge(v, u, mesh.get());
  }
  /*!
   * 			Function to clear edge'face incidency and vertex incidency if
   * present and then to remove the involved edge from the mesh. Resulting
   *       isolated vertices are removed.
   * \param	v	The first vertex of the involved edge.
   * \param	u	The second vertex of the involved edge.
   * \param	mesh	The involving mesh.
   */
  static void
  remove_edge(vertex_descriptor v, vertex_descriptor u, ref_mesh mesh)
  {
    remove_edge(v, u, &mesh);
  }
  /*!
   * 			Function to clear face'edge incidency and vertex incidency if
   * present and then to remove the involved face from the mesh. Resulting
   *       isolated edges are removed.
   * \param	face	The involving face.
   * \param	mesh	The involving mesh.
   */
  static void remove_face(face_descriptor face, ptr_mesh mesh)
  {
    if(face == null_face())
      throw std::invalid_argument(
          "AIFTopologyHelpers::remove_face(f, m) -> face is null, cannot "
          "suppress it from mesh.");
    // correction towards CGAL-like implementation [done on 2017-05-03]
    // removes edges from the graph if they were already border edges.
    // If this creates isolated vertices remove them as well (at the the
    // remove_edge level).
    auto edgesRange = face->GetIncidentEdges();
    std::list< edge_descriptor > edgesToRemove;
    std::set< vertex_descriptor > vertexWhichNeedsToClearPrecomputation;
    auto itE = edgesRange.begin();
    for(; itE != edgesRange.end(); ++itE)
    {
      vertexWhichNeedsToClearPrecomputation.insert((*itE)->get_first_vertex());
      vertexWhichNeedsToClearPrecomputation.insert((*itE)->get_second_vertex());
      if((*itE)->GetDegree() == 1)
      {
        edgesToRemove.push_back(*itE);
      }
    }
    unlink_all_edges(face);
    mesh->EraseIsolatedFace(face);
    std::list< edge_descriptor >::iterator it = edgesToRemove.begin(),
                                           ite = edgesToRemove.end();
    for(; it != ite; ++it)
    {
      unlink_all_vertices(*it);
      remove_edge(*it, mesh);
    }

    std::set< vertex_descriptor >::iterator
        itv = vertexWhichNeedsToClearPrecomputation.begin(),
        itve = vertexWhichNeedsToClearPrecomputation.end();
    for(; itv != itve; ++itv)
    {
      (*itv)->m_Incident_PtrFaces.clear();
      (*itv)->m_Incident_PtrFaces_Computed = false;
    }
  }
  /*!
   * 			Function to clear face'edge incidency and vertex incidency if
   * present and then to remove the involved face from the mesh. Resulting
   *       isolated edges are removed.
   * \param	face	The involving face.
   * \param	mesh	The involving mesh.
   */
  static void remove_face(face_descriptor face, smart_ptr_mesh mesh)
  {
    remove_face(face, mesh.get());
  }
  /*!
   * 			Function to clear face'edge incidency and vertex incidency if
   * present and then to remove the involved face from the mesh. Resulting
   *       isolated edges are removed.
   * \param	face	The involving face.
   * \param	mesh	The involving mesh.
   */
  static void remove_face(face_descriptor face, ref_mesh mesh)
  {
    remove_face(face, &mesh);
  }
  /*!
   * Generic function to remove faces given by an iterator
   * range.
   *
   * \param  first  The iterator pointing on the range beginning.
   * \param  last   The iterator pointing on the range end.
   * \param  mesh   The involving mesh.
   */
  template< typename InputIt, typename MeshType >
  static void remove_faces(InputIt first, InputIt last, MeshType mesh)
  {
    for(; first != last; ++first)
      remove_face(*first, mesh);
  }
  /*!
   * Function to unlink the edge' 2 incident vertices and
   * edge.
   *
   * \param  edge  The involving edge.
   */
  static void unlink_all_vertices(edge_descriptor edge)
  {
    vertex_descriptor v1 = edge->get_first_vertex(),
                      v2 = edge->get_second_vertex();
    if(v1 != null_vertex())
      unlink_vertex_and_edge(v1, edge);
    if(v2 != null_vertex())
      unlink_vertex_and_edge(v2, edge);
  }
  /*!
   * 			Function to unlink the vertex'incident edges and vertex.
   * \param	vertex	The involving vertex.
   */
  static void unlink_all_edges(vertex_descriptor vertex)
  {
    typedef edge_container_in_vertex::const_iterator it_type;
    boost::iterator_range< it_type > edges_range = incident_edges(vertex);
    std::vector< edge_descriptor > edgesToUnlink;

    std::copy(
        edges_range.begin(),
        edges_range.end(),
        std::back_inserter(edgesToUnlink)); // the relations need to be copied
                                            // to keep iterator valid

    for(std::vector< edge_descriptor >::const_iterator it =
            edgesToUnlink.cbegin();
        it != edgesToUnlink.cend();
        ++it)
      unlink_vertex_and_edge(vertex, *it);
  }
  /*!
   * 			Function to unlink the face'incident edges and face.
   * \param	face	The involving face.
   */
  static void unlink_all_edges(face_descriptor face)
  {
    typedef incident_edge_iterator it_type;
    boost::iterator_range< it_type > edges_range = incident_edges(face);
    std::vector< edge_descriptor > edgesToUnlink;

    std::copy(
        edges_range.begin(),
        edges_range.end(),
        std::back_inserter(edgesToUnlink)); // the relations need to be copied
                                            // to keep iterator valid

    for(std::vector< edge_descriptor >::const_iterator it =
            edgesToUnlink.cbegin();
        it != edgesToUnlink.cend();
        ++it)
      unlink_edge_and_face(*it, face);
  }
  /*!
   * 			Function to unlink the edge'incident faces and edge.
   * \param	face	The involving face.
   */
  static void unlink_all_faces(edge_descriptor edge)
  {
    typedef face_container_in_edge::const_iterator it_type;
    boost::iterator_range< it_type > faces_range = incident_faces(edge);
    std::vector< face_descriptor > facesToUnlink;

    std::copy(
        faces_range.begin(),
        faces_range.end(),
        std::back_inserter(facesToUnlink)); // the relations need to be copied
                                            // to keep iterator valid

    for(std::vector< face_descriptor >::const_iterator it =
            facesToUnlink.cbegin();
        it != facesToUnlink.cend();
        ++it)
      unlink_edge_and_face(edge, *it);
  }

  ////////////////////////////////////////////////////////
  /*!
   * Add the vertex to the edge'incident vertices (at specified
   * position) if edge and vertex were not incident. Update caching information
   * for vertex. Usage precondition: vertex != null_vertex() and edge !=
   * null_edge()
   *
   * \param  edge      The involving edge
   * \param  vertex    The involving vertex
   * \param  position  The position of the vertex in the edge (FIRST or
   *                   SECOND)
   *
   * \see  add_vertex_to_edge(edge,  vertex, position)
   */
  static void add_vertex(edge_descriptor edge,
                         vertex_descriptor vertex,
                         enum vertex_pos position)
  {
    if(vertex == null_vertex())
      throw std::invalid_argument(
          "AIFTopologyHelpers::add_edge(v, e) -> vertex is a null vertex.");

    if(edge == null_edge())
      throw std::invalid_argument(
          "AIFTopologyHelpers::add_edge(v, e) -> edge is a null edge.");

    if(!are_incident(edge, vertex))
    {
      add_vertex_to_edge(edge, vertex, position);
    }
  }
  /*!
   * Remove the vertex from the edge'incident vertices if edge and
   * vertex were incident. Update caching information for vertex. Usage
   * precondition: edge != null_edge() and vertex != null_vertex()
   *
   * \param  edge    The edge for which incident vertices are reduced.
   * \param  vertex  The vertex to remove.
   *
   * \see  remove_vertex_from_edge(edge,  vertex)
   */
  static void remove_vertex(edge_descriptor edge, vertex_descriptor vertex)
  {
    if(edge == null_edge())
      throw std::invalid_argument(
          "AIFTopologyHelpers::remove_vertex(e, v) -> edge is a null edge.");

    if(vertex == null_vertex())
      throw std::invalid_argument("AIFTopologyHelpers::remove_vertex(e, v) -> "
                                  "vertex is a null vertex.");

    if(are_incident(edge, vertex))
    {
      remove_vertex_from_edge(edge, vertex);
    }
  }
  /*!
   * Add the edge to the vertex'incident edges if vertex and edge were
   * not incident. Update caching information for vertex' incident edges. Usage
   * precondition: vertex != null_vertex() and edge != null_edge()
   *
   * \param  vertex  The vertex for which incident edges are augmented.
   * \param  edge    The edge to add.
   *
   * \see  add_edge_to_vertex(vertex, edge)
   */
  static void add_edge(vertex_descriptor vertex, edge_descriptor edge)
  {
    if(vertex == null_vertex())
      throw std::invalid_argument(
          "AIFTopologyHelpers::add_edge(v, e) -> vertex is a null vertex.");

    if(edge == null_edge())
      throw std::invalid_argument(
          "AIFTopologyHelpers::add_edge(v, e) -> edge is a null edge.");

    if(!are_incident(vertex, edge))
    {
      add_edge_to_vertex(vertex, edge);
    }
  }
  /*!
   * 			Remove the edge from the vertex'incident edges if vertex and edge
   * were incident. Update caching information for face' incident vertices.
   * Usage precondition: vertex != null_vertex() and edge != null_edge()
   * \param	edge	The edge for which incident faces are reduced.
   * \param	face	The face to remove.
   * \see remove_edge_from_vertex(vertex, edge)
   */
  static void remove_edge(vertex_descriptor vertex, edge_descriptor edge)
  {
    if(vertex == null_vertex())
      throw std::invalid_argument(
          "AIFTopologyHelpers::remove_edge(v, e) -> vertex is a null vertex.");
    if(edge == null_edge())
      throw std::invalid_argument(
          "AIFTopologyHelpers::remove_edge(v, e) -> edge is a null edge.");

    if(are_incident(vertex, edge))
    {
      remove_edge_from_vertex(vertex, edge);
    }
  }
  /*!
   * 			Add the face to the edge'incident faces if edge and face were not
   * incident. Update caching information for face' incident vertices. Usage
   * precondition: edge != null_edge() and face != null_face()
   * \param	edge	The edge for which incident faces are augmented.
   * \param	face	The face to add.
   * \see add_face_to_edge(edge, face)
   */
  static void add_face(edge_descriptor edge, face_descriptor face)
  {
    if(edge == null_edge())
      throw std::invalid_argument(
          "AIFTopologyHelpers::add_face(e, f) -> edge is a null edge.");
    if(face == null_face())
      throw std::invalid_argument(
          "AIFTopologyHelpers::add_face(e, f) -> face is a null face.");

    if(!are_incident(edge, face))
    {
      add_face_to_edge(edge, face);
    }
  }
  /*!
   * 			Remove the face from the edge'incident faces if edge and face were
   * incident. Update caching information for face' incident vertices. Usage
   * precondition: edge != null_edge() and face != null_face()
   * \param	edge	The edge for which incident faces are reduced.
   * \param	face	The face to remove.
   * \see remove_face_from_edge(edge, face)
   */
  static void remove_face(edge_descriptor edge, face_descriptor face)
  {
    if(edge == null_edge())
      throw std::invalid_argument(
          "AIFTopologyHelpers::remove_face(e, f) -> edge is a null edge.");
    if(face == null_face())
      throw std::invalid_argument(
          "AIFTopologyHelpers::remove_face(e, f) -> face is a null face.");

    if(are_incident(edge, face))
    {
      remove_face_from_edge(edge, face);
    }
  }
  /*!
   * 			Add the edge to the face'incident edges if face and edge were not
   * incident. Update caching information for face' incident vertices. Usage
   * precondition: face != null_face() and edge != null_edge()
   * \param	face	The face for which incident edges are augmented.
   * \param	edge	The edge to add.
   * \see add_edge_to_face(face, edge)
   */
  static void add_edge(face_descriptor face, edge_descriptor edge)
  {
    if(edge == null_edge())
      throw std::invalid_argument(
          "AIFTopologyHelpers::add_edge(f, e) -> edge is a null edge.");
    if(face == null_face())
      throw std::invalid_argument(
          "AIFTopologyHelpers::add_edge(f, e) -> face is a null face.");

    if(!are_incident(face, edge))
    {
      add_edge_to_face(face, edge);
    }
  }
  /*!
   * 			Remove the edge to the face'incident edges if face and edge were
   * incident. Update caching information for face' incident vertices. Usage
   * precondition: face != null_face() and edge != null_edge()
   * \param	face	The face for which incident edges are reduced.
   * \param	edge	The edge to remove.
   * \see remove_edge_from_face(face, edge)
   */
  static void remove_edge(face_descriptor face, edge_descriptor edge)
  {
    if(edge == null_edge())
      throw std::invalid_argument(
          "AIFTopologyHelpers::add_edge(f, e) -> edge is a null edge.");
    if(face == null_face())
      throw std::invalid_argument(
          "AIFTopologyHelpers::add_edge(f, e) -> face is a null face.");

    if(are_incident(face, edge))
    {
      remove_edge_from_face(face, edge);
    }
  }

public:
  // One-ring stuff
  /*!
   *           Return the full one-ring of a vertex,
   *           not restricted to edges having their 2
   *           vertices adjacent to v.
   */
  static std::set< edge_descriptor >
  get_unordered_one_ring_edges(vertex_descriptor v)
  {
    std::set< edge_descriptor > edgesToRemove(incident_edges(v).begin(),
                                              incident_edges(v).end());
    std::set< edge_descriptor > edges;
    std::set< edge_descriptor > oldEdges;
    bool first = true;
    ////////////////////////////////////////////////////////////////////////////
    face_container_in_vertex faces = incident_faces_container(v);
    face_container_in_vertex::iterator nextFace = faces.begin();
    while(!faces.empty())
    {
      auto eRange = incident_edges(*nextFace);
      edges.insert(eRange.begin(), eRange.end());

      if(first)
      {
        first = false;
        ////////////////////////////////////////////////////////////////////////////
        oldEdges = edges;
      }
      else
      {
        std::set< edge_descriptor > edgesTmp(eRange.begin(), eRange.end());
        std::set< edge_descriptor > edgesInter;
        std::set_intersection(edgesTmp.begin(),
                              edgesTmp.end(),
                              oldEdges.begin(),
                              oldEdges.end(),
                              std::inserter(edgesInter, edgesInter.end()));
        oldEdges.clear();
        std::set_difference(edges.begin(),
                            edges.end(),
                            edgesInter.begin(),
                            edgesInter.end(),
                            std::inserter(oldEdges, oldEdges.end()));
        ////////////////////////////////////////////////////////////////////////////
        edges.swap(oldEdges); // edges = oldEdges;
      }
      ////////////////////////////////////////////////////////////////////////////
      faces.erase(nextFace); // each processed face is removed
      ////////////////////////////////////////////////////////////////////////////
      nextFace = faces.begin();
    }
    ////////////////////////////////////////////////////////////////////////////
    std::set< edge_descriptor > result;
    std::set< edge_descriptor >::iterator it(edges.begin()), ite(edges.end());
    for(; it != ite; ++it)
    {
      if(edgesToRemove.find(*it) ==
         edgesToRemove.end()) // *it is not an incident edge
      {
        result.insert(*it);
      }
    }
    return result;
  }

  /**
   * return true when edge e is an element of container_of_edge_one_ring whose
   * two vertices are incident to only one (resp. possible many) other edge of
   * container_of_edge_one_ring if allow_multiple_incident_edges is false (resp.
   * true).
   */
  template< typename T >
  static bool
  is_e_one_ring_connected(edge_descriptor e,
                          const T &container_of_edge_one_ring,
                          bool allow_multiple_incident_edges = false, 
                          bool ensure_e_belongs_to_container = true)
  {
    // First ensure that e is an element of container_of_edge_one_ring
    if(ensure_e_belongs_to_container && 
      (std::find(container_of_edge_one_ring.begin(),
                 container_of_edge_one_ring.end(),
                 e) == container_of_edge_one_ring.end()))
      return false;
    // Then check if e is one-ring connected to its one ring
    bool another_edge_incident_to_e_v1 = false,
         another_edge_incident_to_e_v2 = false;
    for(auto element : container_of_edge_one_ring)
    {
      if((element != e) && are_incident(element, e->get_first_vertex()))
      {
        if(!allow_multiple_incident_edges)
          if(another_edge_incident_to_e_v1)
            return false; // found 2 one-ring edges incident to v1
        another_edge_incident_to_e_v1 = true;
        if(allow_multiple_incident_edges)
          if(another_edge_incident_to_e_v2)
            break;
      }
      else if((element != e) && are_incident(element, e->get_second_vertex()))
      {
        if(!allow_multiple_incident_edges)
          if(another_edge_incident_to_e_v2)
            return false; // found 2 one-ring edges incident to v2
        another_edge_incident_to_e_v2 = true;
        if(allow_multiple_incident_edges)
          if(another_edge_incident_to_e_v1)
            break;
      }
    }

    return another_edge_incident_to_e_v1 && another_edge_incident_to_e_v2;
  }

  /**
   * return true if all adjacent edges to e that are also elements of
   * container_of_edge_one_ring are one ring connected (using function
   * is_e_one_ring_connected).
   */
  template< typename T >
  static bool
  are_adjacent_edges_one_ring_connected(edge_descriptor e,
                                        const T &container_of_edge_one_ring)
  {
    auto edge_range = adjacent_edges(e);
    auto iter = edge_range.begin();
    for(; iter != edge_range.end(); ++iter)
    {
      if((std::find(container_of_edge_one_ring.begin(),
                    container_of_edge_one_ring.end(),
                    *iter) != container_of_edge_one_ring.end()) &&
         !is_e_one_ring_connected(*iter, container_of_edge_one_ring))
        return false;
    }
    return true;
  }
  /**
   * return the first edge of container_of_edge_one_ring
   * which satisfies the are_adjacent_edges_one_ring_connected(e,
   * container_of_edge_one_ring) predicate.
   */
  template< typename T >
  static edge_descriptor
  get_first_edge_whose_adjacent_edges_are_one_ring_connected(
      const T &container_of_edge_one_ring)
  {
    for(auto element : container_of_edge_one_ring)
      if(are_adjacent_edges_one_ring_connected(element,
                                               container_of_edge_one_ring))
        return element;
    // else return the first one_ring_connected edge (without multiple incident
    // edges)
    for(auto element : container_of_edge_one_ring)
      if(is_e_one_ring_connected(element, container_of_edge_one_ring, false))
        return element;
    // else return the first one_ring_connected edge (with multiple incident
    // edges)
    for(auto element : container_of_edge_one_ring)
      if(is_e_one_ring_connected(element, container_of_edge_one_ring, true))
        return element;
    // else return the first element of the container
    return *(container_of_edge_one_ring.begin());
  }
  /*!
   *           Return the full one-ring of a vertex,
   *           not restricted to adjacent edges.
   *           Always ordered in clockwise order when 
   *           faces are oriented into counter-clockwise.
   */
  static edge_container_in_vertex
  get_ordered_one_ring_edges(vertex_descriptor v)
  {
    std::set< edge_descriptor > notSortedOneRing =
        get_unordered_one_ring_edges(v);
    std::set< edge_descriptor > notSortedOneRingSave(notSortedOneRing);
    edge_container_in_vertex result;

    if(notSortedOneRing.size() == 0)
      return result;

    result.reserve(notSortedOneRing.size());

    edge_descriptor nextEdge =
        get_first_edge_whose_adjacent_edges_are_one_ring_connected(
            notSortedOneRing); // this call is needed to ensure that the reverse
                               // step works properly Indeed, this needs that
                               // the first one-ring edge is "classical"
                               // one-ring edge
    while(!notSortedOneRing.empty())
    {
      result.push_back(nextEdge);
      notSortedOneRing.erase(nextEdge);
      if(!notSortedOneRing.empty())
      {
        if(notSortedOneRing.size() == 1)
        {
          nextEdge = *notSortedOneRing.begin();
        }
        else
        {
          std::set< edge_descriptor >::iterator it(notSortedOneRing.begin()),
              ite(notSortedOneRing.end());
          for(; it != ite; ++it)
            if(are_adjacent(*it, nextEdge))
            {
              // next we check for a non-manifold configuration (e.g. a
              // complex edge or a dangling edge here)
              // 1) we search for another edge satisfying the condition
              // are_adjacent(*it, nextEdge) 2) if such another edge exists,
              // then
              //    => among all these edge(>=2) we select the first which has
              //    its second vertex not incident
              //       to any of one-ring edges (remaining notSortedOneRing and
              //       already used ones)
              //    => If this fails, look for an edge adjacent to last edge 
              //    with lowest degree for common vertex
              // 3) if there are only edge(>=2) with second vertex incident to
              // one of remaining edge, then
              //    select the first
              // Remark: for 2 edges in the same situation, selecting the first
              // is because without geometric information we cannot propose
              // something else.
              ///////////////////////////////////////////////////////////////
              edge_descriptor nextEdgeTmp = nextEdge;
              nextEdge = *it;
              // 1)
              auto it2 = it;
              ++it2;

              bool neverEnter = true, atLeastOneAnotherAdjacent = false;
              for(; it2 != ite; ++it2)
              {
                if (are_adjacent(*it2, nextEdgeTmp))
                {
                  atLeastOneAnotherAdjacent = true;
                  if (are_adjacent(*it2, nextEdge))
                  {
                    neverEnter = false;
                    // 2)
                    vertex_descriptor common_v = common_vertex(nextEdgeTmp, *it2);
                    vertex_descriptor other_v = opposite_vertex(*it2, common_v);
                    auto it3 = it;
                    ++it3;
                    for (; it3 != ite; ++it3)
                      if ((it3 != it2) && are_incident(*it3, other_v))
                      {
                        break; /// \todo  Should not break in case of a dangling polyline
                      }
                    if ((it3 == ite) && (!is_e_one_ring_connected(*it2, notSortedOneRingSave, true, true) ||
                      (is_e_one_ring_connected(nextEdge, notSortedOneRingSave, true, true) && 
					  is_e_one_ring_connected(*it2, result, true, false))))
                    {
                      nextEdge = *it2;
                      break;
                    }
                  }
                }
              }
              ///////////////////////////////////////////////////////////////
              if (atLeastOneAnotherAdjacent && neverEnter)
              {
                it2 = it;
                ++it2;

                for (; it2 != ite; ++it2)
                {
                  if (are_adjacent(*it2, nextEdgeTmp) && degree(common_vertex(nextEdgeTmp, nextEdge))>degree(common_vertex(nextEdgeTmp, *it2)))
                  {
                    nextEdge = *it2;
                    break;
                  }
                }
              }
              ///////////////////////////////////////////////////////////////
              break;
            }

          if(it == ite)
          {
            // try to go back at the beginning in the other direction
            std::set< edge_descriptor > notUsedEdgesCopy(notSortedOneRing);

            edge_descriptor nextTmpEdge =
                *result.begin(); // get the first used edge
            it = notUsedEdgesCopy.begin();
            ite = notUsedEdgesCopy.end();
            while(it != ite)
            {
              bool found = false;
              for(; it != ite; ++it)
                if(are_adjacent(*it, nextTmpEdge))
                {
                  nextTmpEdge = *it;
                  notUsedEdgesCopy.erase(*it);
                  found = true;
                  break;
                }

              if(!found)
                break;
              it = notUsedEdgesCopy.begin();
              ite = notUsedEdgesCopy.end();
            }

            // There are 2 cases to take into account: 1) an hole, while the
            // vertex this is still 2-manifold (border case); 2) several
            // components and the vertex this is a cut vertex
            if(is_cut_vertex(v)                // non-manifold
               || notUsedEdgesCopy.size() > 0) // hole within the one-ring
            { // non-manifold case: not all cases managed => to do later
              if (notUsedEdgesCopy.size() == 0)
              {                
                if (notSortedOneRing.size() == 2)
                {
                  if (are_adjacent(*notSortedOneRing.begin(), result.front()))
                  {
                    auto it_tmp = notSortedOneRing.begin();
                    ++it_tmp;
                    nextEdge = *it_tmp;
                  }
                  else
                    nextEdge = *notSortedOneRing.begin();
                }
                else
                  nextEdge = *notSortedOneRing.begin();
              }
              else if(notUsedEdgesCopy.size() <= 2)
                nextEdge = *notUsedEdgesCopy.begin();
              else
              {
                nextTmpEdge = *notUsedEdgesCopy.begin();
                it = notUsedEdgesCopy.begin();
                ite = notUsedEdgesCopy.end();
                while(it != ite)
                {
                  bool found = false;
                  for(; it != ite; ++it)
                    if(are_adjacent(*it, nextTmpEdge))
                    {
                      nextTmpEdge = *it;
                      notUsedEdgesCopy.erase(*it);
                      found = true;
                      break;
                    }
                  if(!found)
                    break;
                  it = notUsedEdgesCopy.begin();
                  ite = notUsedEdgesCopy.end();
                }

                nextEdge = nextTmpEdge;
              }
            }
            else
              nextEdge = nextTmpEdge;
          }
        }
      }
    }
    // Try to get a one-ring ordered into clockwise order when faces are
    // oriented into counter-clockwise order The code bellow is needed when the
    // first not ordered one-ring edge is followed by a hole in the clockwise
    // order
    size_t nbE = result.size();
    if(nbE > 2)
    {
      size_t i;
      face_container_in_vertex faces;
      for(i = 0; i < nbE; ++i)
      {
        if(are_adjacent(result[i], result[(i + 1) % nbE]))
        {
          auto fRange = incident_faces(result[i]); // face iterator range
          if(!fRange.empty())
          {
            faces = face_container_in_vertex(fRange.begin(), fRange.end());
            break;
          }
        }
      }
      if(i < nbE)
      {
        // try to find a face that has "this" vertex as incident vertex
        vertex_container_in_vertex resV;
        size_t nbF = faces.size();
        size_t k;
        for(k = 0; k < nbF; ++k)
        {
          if(are_incident(v, faces[k]))
          {
            auto vRange =
                incident_vertices(faces[k]); // vertices iterator range
            resV = vertex_container_in_vertex(vRange.begin(), vRange.end());
			break;
          }
        }
        if(k == nbF)
        {
          auto vRange = incident_vertices(faces[0]); // vertices iterator range
          resV = vertex_container_in_vertex(vRange.begin(), vRange.end());
        }
        //////////////////////////////////////////////////////////////////////////
        // the order generally given for face in mesh files is the
        // counterclockwise order
        bool v1BeforeV2InFace = false;
        //////////////////////////////////////////////////////////////////////////
        auto it = resV.begin(); // iterator over a container of pointers
        auto ite = resV.end();
        vertex_descriptor first_v, second_v;
        if(are_adjacent(result[i], result[(i + 1) % nbE]))
          second_v = common_vertex(result[i], result[(i + 1) % nbE]);
        else
          second_v = result[i]->get_second_vertex();
        first_v = opposite_vertex(result[i], second_v);
        while(it != ite)
        {
          if(**it == *(first_v))
          {
            ++it;
            if(it == ite)
              it = resV.begin();
            if(**it == *(second_v))
              v1BeforeV2InFace =
                  true; // if the following is V2 => true, else => false

            break;
          }
          else if(**it == *(second_v))
          {
            ++it;
            if (it == ite)
              it = resV.begin();
            if(**it != *(first_v))
              v1BeforeV2InFace = true; // if the following is different from V1
                                       // => true, else => false

            break;
          }
          ++it;
        }
        if(it != ite)
        { // ensure that face orientation (given by its vertex indices order) is
          // opposite to one-ring orientation
          if(are_incident(first_v, result[(i + 1) % nbE]))
          {
            if(!v1BeforeV2InFace) // in that case currently the one-ring is
                                  // clockwise oriented => change it to
                                  // counterclockwise!
            {
              std::reverse(result.begin(), result.end());
            }
          }
          else if(are_incident(second_v, result[(i + 1) % nbE]))
          {
            if(v1BeforeV2InFace) // in that case currently the one-ring is
                                 // counter clockwise oriented => change it to
                                 // clockwise!
            {
              std::reverse(result.begin(), result.end());
            }
          }
        }
      }
    }
    return result;
  }

  /*!
   *           Return the full one-ring of a vertex,
   *           not restricted to adjacent vertices.
   *           Always ordered in clockwise order when 
   *           faces are oriented into counter-clockwise.
   */
  static vertex_container_in_vertex
  get_ordered_one_ring_vertices(vertex_descriptor v)
  {
    auto &m_One_Ring_Vertices = v->m_One_Ring_Vertices;

    if(v->m_Is_One_Ring_Vertices_Computed)
    {
      // DBG std::cout << "Use one-ring-vertices cache" << std::endl;
      return m_One_Ring_Vertices; // already in cache
    }

    ////////////////////////////////////////////////////////////////////////////
    // Currently compute the ordered one-ring of vertices
    // (not necessarily adjacent vertices if some incident faces are not
    // triangular)
    // DBG std::cout << "No valid one-ring-vertices cache, computing it" <<
    // std::endl;
    m_One_Ring_Vertices.clear();
    edge_container_in_vertex orderedOneRing = get_ordered_one_ring_edges(v);
    vertex_container_in_vertex oneR;
    size_t nbE = orderedOneRing.size();
    if(nbE == 0)
      return m_One_Ring_Vertices;
    oneR.reserve(nbE);
    if(nbE == 1)
    {
      if(orderedOneRing[0]->get_first_vertex() ==
         orderedOneRing[0]->get_second_vertex())
      {
        oneR.push_back(orderedOneRing[0]->get_first_vertex());
      }
      else
      {
        oneR.push_back(orderedOneRing[0]->get_first_vertex());
        oneR.push_back(orderedOneRing[0]->get_second_vertex());
      }
    }
    else
    {
      vertex_descriptor lastVertex;

      if(are_incident(orderedOneRing[0]->get_second_vertex(),
                      orderedOneRing[1]))
      {
        oneR.push_back(orderedOneRing[0]->get_first_vertex());
        if(orderedOneRing[0]->get_first_vertex() !=
           orderedOneRing[0]->get_second_vertex())
          oneR.push_back(orderedOneRing[0]->get_second_vertex());
        lastVertex = orderedOneRing[0]->get_second_vertex();
      }
      else if(are_incident(orderedOneRing[0]->get_first_vertex(),
                           orderedOneRing[1]))
      {
        oneR.push_back(orderedOneRing[0]->get_second_vertex());
        if(orderedOneRing[0]->get_first_vertex() !=
           orderedOneRing[0]->get_second_vertex())
          oneR.push_back(orderedOneRing[0]->get_first_vertex());
        lastVertex = orderedOneRing[0]->get_first_vertex();
      }
      else
      {
        // the 2 one-ring edges are not connected
        oneR.push_back(orderedOneRing[0]->get_first_vertex());
        if(orderedOneRing[0]->get_first_vertex() !=
           orderedOneRing[0]->get_second_vertex())
          oneR.push_back(orderedOneRing[0]->get_second_vertex());
        lastVertex = orderedOneRing[0]->get_second_vertex();
      }

      for(size_t i = 1; i < nbE; ++i)
      {
        if(are_incident(lastVertex, orderedOneRing[i]))
        {
          lastVertex = opposite_vertex(orderedOneRing[i], lastVertex);
        }
        else
        {
          if(i == nbE - 1)
          {
            if(are_incident(oneR[0], orderedOneRing[i]))
            {
              oneR.push_back(opposite_vertex(orderedOneRing[i], oneR[0]));
            }
            else if(are_incident(oneR[1], orderedOneRing[i]))
            {
              oneR.push_back(opposite_vertex(orderedOneRing[i], oneR[1]));
            }
            else
            {
              oneR.push_back(orderedOneRing[i]->get_first_vertex());
              if(orderedOneRing[i]->get_first_vertex() !=
                 orderedOneRing[i]->get_second_vertex())
                oneR.push_back(orderedOneRing[i]->get_second_vertex());
            }

            break;
          }
          else
          {
            if(are_incident(orderedOneRing[i]->get_second_vertex(),
                            orderedOneRing[i + 1]))
            {
              if(std::find(oneR.begin(),
                           oneR.end(),
                           orderedOneRing[i]->get_second_vertex()) !=
                 oneR.end())
              {
                if(orderedOneRing[i]->get_first_vertex() !=
                   orderedOneRing[i]->get_second_vertex())
                {
                  oneR.push_back(orderedOneRing[i]->get_second_vertex());
                }
                lastVertex = orderedOneRing[i]->get_first_vertex();
              }
              else
              {
                if(orderedOneRing[i]->get_first_vertex() !=
                   orderedOneRing[i]->get_second_vertex())
                {
                  oneR.push_back(orderedOneRing[i]->get_first_vertex());
                }
                lastVertex = orderedOneRing[i]->get_second_vertex();
              }
            }
            else if(are_incident(orderedOneRing[i]->get_first_vertex(),
                                 orderedOneRing[i + 1]))
            {
              if(std::find(oneR.begin(),
                           oneR.end(),
                           orderedOneRing[i]->get_first_vertex()) != oneR.end())
              {
                if(orderedOneRing[i]->get_first_vertex() !=
                   orderedOneRing[i]->get_second_vertex())
                {
                  oneR.push_back(orderedOneRing[i]->get_first_vertex());
                }
                lastVertex = orderedOneRing[i]->get_second_vertex();
              }
              else
              {
                if(orderedOneRing[i]->get_first_vertex() !=
                   orderedOneRing[i]->get_second_vertex())
                {
                  oneR.push_back(orderedOneRing[i]->get_second_vertex());
                }
                lastVertex = orderedOneRing[i]->get_first_vertex();
              }
            }
            else
            {
              // the 2 one-ring edges are not connected
              if(orderedOneRing[i]->get_first_vertex() !=
                 orderedOneRing[i]->get_second_vertex())
              {
                oneR.push_back(orderedOneRing[i]->get_first_vertex());
              }
              lastVertex = orderedOneRing[i]->get_second_vertex();
            }
          }
        }

        if(std::find(oneR.begin(), oneR.end(), lastVertex) != oneR.end())
        {
          if ((i + 1 == nbE) && !is_e_one_ring_connected(
                                orderedOneRing[0], orderedOneRing, true))
          {
            oneR.push_back(lastVertex);
            continue;
          }
          else if((i + 2 == nbE) && !is_e_one_ring_connected(
                                   orderedOneRing[i + 1], orderedOneRing, true))
          {
            oneR.push_back(lastVertex);
            continue;
          }
          break;
        }

        oneR.push_back(lastVertex);
      }
    }

    m_One_Ring_Vertices = vertex_container_in_vertex(oneR.begin(), oneR.end());
    v->m_Is_One_Ring_Vertices_Computed = true;

    return m_One_Ring_Vertices;
  }

  /*!
   * 			Return the one-ring of adjacent vertices.
   *           Always ordered in clockwise order when 
   *           faces are oriented into counter-clockwise.
   */
  static vertex_container_in_vertex
  get_ordered_one_ring_of_adjacent_vertices(vertex_descriptor v)
  {
    vertex_container_in_vertex resFullRing = get_ordered_one_ring_vertices(
        v); // ring of vertices not necessarily adjacent
    vertex_container_in_vertex oneR;
    size_t nbE = resFullRing.size();
    if(nbE == 0)
      return oneR;

    oneR.reserve(nbE);
    ///////////////////////////////////////////////////////////////////////////
    auto it(resFullRing.begin());
    auto ite(resFullRing.end());
    while(it != ite)
    {
      if(are_adjacent(*it, v))
        oneR.push_back(*it);

      ++it;
    }
    /////////////////////////////////////////////////////////////////////////
    return vertex_container_in_vertex(oneR.begin(), oneR.end());
  }
  /*!
   * 			A 2 manifold vertex one ring
   * \param	v	The involving vertex
   * \return	A boolean telling if the one-ring of v is 2-manifold, considering
   * only the topology (not the geometry).
   */
  static bool is_one_ring_2_manifold(vertex_descriptor v)
  {
    if(!is_2_manifold_vertex(v))
      return false;

    std::set< edge_descriptor > starEdges = get_unordered_one_ring_edges(v);
    auto its = starEdges.begin();
    auto itse = starEdges.end();
    for(; its != itse; ++its)
    {
      if(is_complex_edge(*its))
        return false;
    }

    return true;
  }
  // End of one-ring stuff
  /////////////////////////////////
  /*!
  * 			Get the edge of the face'incident edges after the prev_edge
  * \param	face	The face for which incident edges are retrieved.
  * \param	prev_edge	The edge located before the returned edge.
  */
  static edge_descriptor get_edge_of_face_after_edge(face_descriptor face, 
                                                     edge_descriptor prev_edge)
  {
    std::vector< edge_descriptor >::iterator
      it = face->m_Incident_PtrEdges.begin(),
      ite = face->m_Incident_PtrEdges.end();
    for (; it != ite; ++it)
    {
      if ((((*it)->get_first_vertex() == prev_edge->get_first_vertex()) &&
        ((*it)->get_second_vertex() == prev_edge->get_second_vertex())) ||
        (((*it)->get_first_vertex() == prev_edge->get_second_vertex()) &&
        ((*it)->get_second_vertex() == prev_edge->get_first_vertex())))
      {
        ++it;
        if (it == ite)
          it = face->m_Incident_PtrEdges.begin();
        break;
      }
    }

    if (it == ite)
      throw std::invalid_argument("Helpers::get_edge_of_face_after_edge(f, pe) -> prev_edge does not belong to input face.");

    return *it;
  }
  /*!
  * 			Get the edge of the face'incident edges before the next_edge
  * \param	face	The face for which incident edges are retrieved.
  * \param	next_edge	The edge located after the returned edge.
  */
  static edge_descriptor get_edge_of_face_before_edge(face_descriptor face, 
                                                      edge_descriptor next_edge)
  {
    std::vector< edge_descriptor >::iterator
      it = face->m_Incident_PtrEdges.begin(),
      ite = face->m_Incident_PtrEdges.end();
    for (; it != ite; ++it)
    {
      if ((((*it)->get_first_vertex() == next_edge->get_first_vertex()) && 
        ((*it)->get_second_vertex() == next_edge->get_second_vertex())) ||
        (((*it)->get_first_vertex() == next_edge->get_second_vertex()) && 
        ((*it)->get_second_vertex() == next_edge->get_first_vertex())))
      {
        if (it == face->m_Incident_PtrEdges.begin())
        {
          it = ite;
        }
        --it;
        break;
      }
    }

    if (it == ite)
      throw std::invalid_argument("Helpers::get_edge_of_face_before_edge(f, ne) -> next_edge does not belong to input face.");

    return *it;
  }

  /*!
   * 			Add the edge to the face'incident edges after the prev_edge
   * without any check. Update caching information for face' incident vertices.
   * \param	face	The face for which incident edges are augmented.
   * \param	prev_edge	The edge after which the edge is added.
   * \param	edge	The edge to add.
   * \see add_edge(face, edge)
   */
  static void add_edge_to_face_after_edge(face_descriptor face,
                                          edge_descriptor prev_edge,
                                          edge_descriptor edge)
  {
    std::vector< edge_descriptor >::iterator
        it = face->m_Incident_PtrEdges.begin(),
        ite = face->m_Incident_PtrEdges.end();
    for(; it != ite; ++it)
    {
      if((((*it)->get_first_vertex() == prev_edge->get_first_vertex()) &&
          ((*it)->get_second_vertex() == prev_edge->get_second_vertex())) ||
         (((*it)->get_first_vertex() == prev_edge->get_second_vertex()) &&
          ((*it)->get_second_vertex() == prev_edge->get_first_vertex())))
      {
        ++it;
        break;
      }
    }
    face->m_Incident_PtrEdges.insert(it, edge);

    face->clear_vertex_incidency(); // clear the caching incidency relation with
                                    // vertices

    // invalidate one-ring-vertices cache
    auto vRange = incident_vertices(face);
    for(auto vIt = vRange.begin(); vIt != vRange.end(); ++vIt)
    {
      if((*vIt)->m_Is_One_Ring_Vertices_Computed)
      {
        (*vIt)->m_Is_One_Ring_Vertices_Computed = false;
        (*vIt)->m_One_Ring_Vertices.clear(); // free memory
        (*vIt)->m_Incident_PtrFaces_Computed = false;
        (*vIt)->m_Incident_PtrFaces.clear(); // free memory
      }
    }
  }

private:
  /*!
   * An incidence relation is created between the vertex and the edge
   * by adding the vertex at the specified position in the edge. Update caching
   * information for vertex.
   *
   * \param  edge      The involving edge
   * \param  vertex    The involving vertex
   * \param  position  The position of the vertex in the edge (FIRST or
   *                   SECOND)
   *
   * \warning  Be careful, this function will replace the former vertex
   *           at the position \c  position  in the edge. Thus you HAVE to
   *           unlink these two element before, otherwise topological
   *           incoherences will be introduce.
   *
   * \see  link_vertex_and_edge(edge, vertex, position)
   * \see  add_vertex(edge, vertex, position)
   */
  static void add_vertex_to_edge(
      edge_descriptor edge,
      vertex_descriptor vertex,
      enum vertex_pos position) // without telling it to the vertex
  {
    if(position == vertex_pos::FIRST)
      edge->set_first_vertex(vertex);
    else
      edge->set_second_vertex(vertex);

    vertex->m_Is_One_Ring_Vertices_Computed = false;
    vertex->m_One_Ring_Vertices.clear(); // free memory
  }
  /*!
   * An incidence relation is created between the vertex and the edge
   * by adding the edge in the incident edges container of the vertex. Update
   * caching information for vertices.
   *
   * \param  vertex  The involving vertex
   * \param  edge    The involving edge
   *
   * \see  link_vertex_and_edge(edge, vertex, position)
   * \see  add_edge(vertex, edge)
   */
  static void add_edge_to_vertex(
      vertex_descriptor vertex,
      edge_descriptor
          edge) // without telling it to the edge (and without verifying that
                // the same edge is already incident to "vertex")
  {
    FEVV::Container::insert(vertex->m_Incident_PtrEdges, edge);

    vertex->m_Is_One_Ring_Vertices_Computed = false;
    vertex->m_One_Ring_Vertices.clear(); // free memory
  }

  /*!
   * An incidence relation is removed between the vertex and the edge
   * by removing the vertex at the specified position in the edge. Update
   * caching information for vertex.
   *
   * \param  edge    The involving edge
   * \param  vertex  The involving vertex
   *
   * \see  unlink_vertex_and_edge(edge, vertex)
   * \see  remove_vertex(edge, vertex)
   */
  static void remove_vertex_from_edge(edge_descriptor edge,
                                      vertex_descriptor vertex)
  {
    if(vertex_position(edge, vertex) == vertex_pos::FIRST)
      edge->set_first_vertex(null_vertex());
    else
      edge->set_second_vertex(null_vertex());

    vertex->m_Is_One_Ring_Vertices_Computed = false;
    vertex->m_One_Ring_Vertices.clear(); // free memory
  }
  /*!
   * An incidence relation is removed between the vertex and the edge
   * by removing the edge from the vertex' incident edges. Update caching
   * information for vertex.
   *
   * \param  vertex  The involving vertex
   * \param  edge    The involving edge
   *
   * \see  unlink_vertex_and_edge(edge, vertex)
   * \see  remove_edge(vertex, edge)
   */
  static void remove_edge_from_vertex(vertex_descriptor vertex,
                                      edge_descriptor edge)
  {
    FEVV::Container::erase(vertex->m_Incident_PtrEdges, edge);

    vertex->m_Is_One_Ring_Vertices_Computed = false;
    vertex->m_One_Ring_Vertices.clear(); // free memory
  }
  /*!
   * Add the edge to the face'incident edges without any
   * check. Update caching information for face' incident vertices.
   *
   * \param  face  The face for which incident edges are augmented.
   * \param  edge  The edge to add.
   *
   * \see  add_edge(face, edge)
   */
  static void add_edge_to_face(face_descriptor face, edge_descriptor edge)
  {
    // DBG std::cout << "add edge " << edge << " to face " << face << std::endl;

    FEVV::Container::insert(face->m_Incident_PtrEdges, edge);

    face->clear_vertex_incidency(); // clear the caching incidency relation with
                                    // vertices

    // invalidate one-ring-vertices cache
    auto vRange = incident_vertices(face);
    for(auto vIt = vRange.begin(); vIt != vRange.end(); ++vIt)
    {
      if((*vIt)->m_Is_One_Ring_Vertices_Computed)
      {
        (*vIt)->m_Is_One_Ring_Vertices_Computed = false;
        (*vIt)->m_One_Ring_Vertices.clear(); // free memory
        (*vIt)->m_Incident_PtrFaces_Computed = false;
        (*vIt)->m_Incident_PtrFaces.clear(); // free memory
      }
    }
  }
  /*!
   * Add the face to the edge'incident faces without any
   * check. Update caching information for edge' incident vertices.
   *
   * \param  edge  The edge for which incident faces are augmented.
   * \param  face  The face to add.
   *
   * \see  add_face(edge,  face)
   */
  static void add_face_to_edge(edge_descriptor edge, face_descriptor face)
  {
    // DBG std::cout << "add face " << face << " to edge " << edge << std::endl;

    FEVV::Container::insert(edge->m_incident_PtrFaces, face);

    // invalidate edge-vertices cache
    auto vRange = incident_vertices(edge);
    for(auto vIt = vRange.begin(); vIt != vRange.end(); ++vIt)
    {
      (*vIt)->m_Is_One_Ring_Vertices_Computed = false;
      (*vIt)->m_One_Ring_Vertices.clear(); // free memory
      (*vIt)->m_Incident_PtrFaces_Computed = false;
      (*vIt)->m_Incident_PtrFaces.clear(); // free memory
    }

    face->clear_vertex_incidency(); // clear the caching incidency relation with
                                    // vertices
  }
  /*!
   * Remove the edge to the face'incident edges without any
   * check. Update caching information for face' incident vertices.
   *
   * \param  face  The face for which incident edges are reduced.
   * \param  edge  The edge to remove.
   *
   * \see  remove_edge(face,  edge)
   */
  static void remove_edge_from_face(face_descriptor face, edge_descriptor edge)
  {
    // DBG std::cout << "remove edge " << edge << " from face " << face <<
    // std::endl;

    // invalidate one-ring-vertices cache
    auto vRange = incident_vertices(face);
    for(auto vIt = vRange.begin(); vIt != vRange.end(); ++vIt)
    {
      (*vIt)->m_Is_One_Ring_Vertices_Computed = false;
      (*vIt)->m_One_Ring_Vertices.clear(); // free memory
      (*vIt)->m_Incident_PtrFaces_Computed = false;
      (*vIt)->m_Incident_PtrFaces.clear(); // free memory
    }

    FEVV::Container::erase(face->m_Incident_PtrEdges, edge);

    face->clear_vertex_incidency(); // clear the caching incidency relation with
                                    // vertices
  }
  /*!
   * Remove the face to the edge'incident faces without any
   * check. Update caching information for edge' incident vertices.
   *
   * \param  edge  The edge for which incident faces are reduced.
   * \param  face  The face to remove.
   *
   * \see  remove_face(edge,  face)
   */
  static void remove_face_from_edge(edge_descriptor edge, face_descriptor face)
  {
    // DBG std::cout << "remove face " << face << " from edge " << edge <<
    // std::endl;

    // invalidate edge-vertices cache
    auto vRange = incident_vertices(edge);
    for(auto vIt = vRange.begin(); vIt != vRange.end(); ++vIt)
    {
      if (*vIt != null_vertex())
      {
        (*vIt)->m_Is_One_Ring_Vertices_Computed = false;
        (*vIt)->m_One_Ring_Vertices.clear(); // free memory
        (*vIt)->m_Incident_PtrFaces_Computed = false;
        (*vIt)->m_Incident_PtrFaces.clear(); // free memory
      }
    }

    FEVV::Container::erase(edge->m_incident_PtrFaces, face);

    face->clear_vertex_incidency(); // clear the caching incidency relation with
                                    // vertices
  }
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
public:
  // HALFEDGE STUFF
  class AIFHalfEdge
  {
  public:
    /// fullfill the DefaultConstructible concept
    AIFHalfEdge(void)
    {
      m_source = null_vertex();
      m_target = null_vertex();
      m_face = null_face();
      m_edge = null_edge();
    }
    /// triplet (source vertex, target vertex, face) constructor
    /// An additionnal halfedge validity check can be done.
    AIFHalfEdge(
        vertex_descriptor s,
        vertex_descriptor t,
        face_descriptor f,
        bool do_validity_check = true) // validity check must be discarded
                                       // during some topological modifications
    {
      // check that the halfedge is valid,
      // aka there is an edge between the 2 vertices...
      std::vector< edge_descriptor > edges = common_edges(s, t);
      if(edges.empty())
        throw std::runtime_error(
            "In AIFHalfEdge::AIFHalfEdge(s, t, f): invalid halfedge, no edge "
            "between the two vertices.");

      m_edge = edges[0];

      if(f == null_face())
      { // appends for for opposite operation on a border halfedge
        m_source =
            s; // we do not want a null m_source to ensure the halfedge exists
        m_target =
            t; // we do not want a null m_target to ensure the halfedge exists
        m_face = f;
        return;
      }
      if(do_validity_check)
      {
        // check that edge belongs to the face
        if(std::none_of(edges.cbegin(), edges.cend(), [f](edge_descriptor e) {
             return are_incident(e, f);
           }))
        {
          throw std::runtime_error(
              "In AIFHalfEdge::AIFHalfEdge(s, t, f): invalid halfedge, not "
              "incident to the face.");
        }
      }
      ///////////////////////////////////////////////////////////////
#if 0 // to ensure that the halfedge is inside the face (not needed if you make
      // sure of that one step before)
				auto edgesRange = incident_edges(f) ;
				auto iter = edgesRange.begin() ;
				edge_descriptor e = *iter;
				++iter;
				edge_descriptor enext = *iter;
				while( !( ((e->get_first_vertex() == s) && (e->get_second_vertex()==t)) || ((e->get_first_vertex() == t) && (e->get_second_vertex()==s))) )
				{
					e = *iter;
					++iter;
					if( iter == edgesRange.end() )
						iter = edgesRange.begin();
					enext = *iter;
				}
				if( e->get_first_vertex() == enext->get_first_vertex()
					|| e->get_first_vertex() == enext->get_second_vertex() )
				{
					m_target = e->get_first_vertex();
					m_source = e->get_second_vertex();
				}
				else if(  e->get_second_vertex() == enext->get_first_vertex()
					|| e->get_second_vertex() == enext->get_second_vertex() )
				{
					m_target = e->get_second_vertex();
					m_source = e->get_first_vertex();
				}
				else
					throw std::runtime_error("In AIFHalfEdge::AIFHalfEdge(s, t, f): unkown case.");
#else
      // the halfedge is not always inside the face (except if you make sure of
      // it one step before)
      m_target = t;
      m_source = s;
#endif
      m_face = f;
    }
    /// quadruplet (source vertex, target vertex, face, corresponding edge)
    /// constructor An additionnal halfedge validity check can be done.
    AIFHalfEdge(
        vertex_descriptor s,
        vertex_descriptor t,
        face_descriptor f,
        edge_descriptor e,
        bool do_validity_check = true) // validity check must be discarded
                                       // during some topological modifications
    {
      m_edge = e;

      if(f == null_face())
      { // appends for for opposite operation on a border halfedge
        m_source =
            s; // we do not want a null m_source to ensure the halfedge exists
        m_target =
            t; // we do not want a null m_target to ensure the halfedge exists
        m_face = f;
        return;
      }
      if(do_validity_check)
      {
        std::vector< edge_descriptor > edges = common_edges(s, t);
        if(edges.empty())
          throw std::runtime_error(
              "In AIFHalfEdge::AIFHalfEdge(s, t, f): invalid halfedge, no edge "
              "between the two vertices.");
        // check that edge belongs to the face
        if(std::none_of(edges.cbegin(), edges.cend(), [f](edge_descriptor e) {
             return are_incident(e, f);
           }))
        {
          throw std::runtime_error(
              "In AIFHalfEdge::AIFHalfEdge(s, t, f): invalid halfedge, not "
              "incident to the face.");
        }
      }
      ///////////////////////////////////////////////////////////////
#if 0 // to ensure that the halfedge is inside the face (not needed if you make
      // sure of that one step before)
				auto edgesRange = incident_edges(f);
				auto iter = edgesRange.begin();
				edge_descriptor edge = *iter;
				++iter;
				edge_descriptor enext = *iter;
				while (!(((edge->get_first_vertex() == s) && (edge->get_second_vertex() == t)) || ((edge->get_first_vertex() == t) && (edge->get_second_vertex() == s))))
				{
					edge = *iter;
					++iter;
					if (iter == edgesRange.end())
						iter = edgesRange.begin();
					enext = *iter;
				}
				if (edge->get_first_vertex() == enext->get_first_vertex()
					|| edge->get_first_vertex() == enext->get_second_vertex())
				{
					m_target = edge->get_first_vertex();
					m_source = edge->get_second_vertex();
				}
				else if (edge->get_second_vertex() == enext->get_first_vertex()
					|| edge->get_second_vertex() == enext->get_second_vertex())
				{
					m_target = edge->get_second_vertex();
					m_source = edge->get_first_vertex();
				}
				else
					throw std::runtime_error("In AIFHalfEdge::AIFHalfEdge(s, t, f): unkown case.");
#else
      // the halfedge is not always inside the face (except if you make sure of
      // it one step before)
      m_target = t;
      m_source = s;
#endif
      m_face = f;
    }

    /// fulfill the EqualityComparable concept
    bool operator==(const AIFHalfEdge &other) const
    {
      // cannot use the incident face, because of some CGAL topological
      // algorithms (the ending condition may be based on an halfedge whose
      // incident face has changed)
      if(m_edge != other.m_edge)
        return false;

      return (
          (m_source == other.m_source) ||
          (m_target == other.m_target)); // sometimes the source is the same,
                                         // sometimes the target is the same
    }

    // SGI (STL) allows != operator based on == operator:
    // https://www.sgi.com/tech/stl/EqualityComparable.html cppreference.com:
    // "For operator!=, the type shall be EqualityComparable."
    bool operator!=(const AIFHalfEdge &other) const
    {
      return !(*this == other);
    }

    /// fulfill the LessThanComparable concept
    bool operator<(const AIFHalfEdge &other) const
    {
      if(m_edge != other.m_edge)
        return m_edge < other.m_edge;
      else
        return m_face < other.m_face;
    }

    /// next of current halfedge
    AIFHalfEdge next(const AIFMesh &m)
    {
      // parse edges around the face
      // look for the edge whose source is current halfedge target
      bool found = false;
      bool isBorderHalfedge = (m_face == null_face());
      if(degree(m_edge) == 0)
      {
        auto edgesRange = get_target()->GetIncidentEdges();
        auto itE = edgesRange.begin();
        for(; itE != edgesRange.end(); ++itE)
        {
          if(are_incident(m_face, *itE))
          {
            return AIFHalfEdge(get_target(),
                               opposite_vertex(*itE, get_target()),
                               m_face,
                               *itE,
                               false);
          }
        }
      }
      edge_descriptor edge /*, edgeAmb*/;
      vertex_type::ptr ptrVn;
      if(isBorderHalfedge)
      {
        // if (is_isolated_vertex(get_target())) // when a modification is in
        // progress (e.g. during the collapse of a border edge) [this approach
        // does not work] 	set_target(opposite_vertex(m_edge, get_source()));

        // Check if target vertex is a 2-manifold vertex
        if(is_2_manifold_vertex(get_target()))
        { // in that case, we know there is only one border edge among
          // all edges incident to m_target
          auto edgesRange = get_target()->GetIncidentEdges();
          auto itE = edgesRange.begin();
          for(; itE != edgesRange.end(); ++itE)
          {
            if(!is_surface_border_edge(*itE))
              continue;

            auto verticesPair = (*itE)->GetIncidentVertices();
            vertex_type::ptr ptrV1 = verticesPair[0];
            vertex_type::ptr ptrV2 = verticesPair[1];

            // we don't know in which order are the vertices in the AIF edge
            if(ptrV1 == get_target() && ptrV2 != get_source())
            {
              found = true;
              ptrVn = ptrV2;
              edge = *itE;
              break;
            }

            if(ptrV2 == get_target() && ptrV1 != get_source())
            {
              found = true;
              ptrVn = ptrV1;
              edge = *itE;
              break;
            }
          }

          if(!found)
          { // particular case due to CGAL low-level topological operations
            itE = edgesRange.begin();
            for(; itE != edgesRange.end(); ++itE)
            {
              auto verticesPair = (*itE)->GetIncidentVertices();
              vertex_type::ptr ptrV1 = verticesPair[0];
              vertex_type::ptr ptrV2 = verticesPair[1];

              // we don't know in which order are the vertices in the AIF edge
              if(*itE == m_edge)
              {
                ++itE; // we take the edge following the edge
                if(itE == edgesRange.end())
                  itE = edgesRange.begin();
                found = true;
                edge = *itE;
                ptrV1 = common_vertex(edge, m_edge);

                return AIFHalfEdge(ptrV1,
                                   opposite_vertex(edge, ptrV1),
                                   incident_faces(edge)[0],
                                   edge,
                                   false);
              }
              else if(ptrV1 == ptrV2)
              {
                found = true;
                edge = *itE;
                return AIFHalfEdge(get_target(),
                                   opposite_vertex(edge, ptrV1),
                                   incident_faces(edge)[0],
                                   edge,
                                   false);
              }
            }
          }
        }
        else if(is_2_manifold_vertex(get_source()))
        { // we try to get the next halfedge by going though
          // the border in the other direction using prev
          AIFHalfEdge h(get_source(), get_target(), m_face, false);
          do
          {
            h = h.prev(m);
          } while(get_target() != h.get_source());
          ptrVn = h.get_target();
          edge = h.m_edge;
          found = true;
        }
        else
          throw std::runtime_error(
              "Error in AIFHalfEdge::next(): next halfedge cannot be found in "
              "a non-manifold context.");
      }
      else
      {
        auto edgesRange = m_face->GetIncidentEdges();
        if((edgesRange.size() == 0) && !are_incident(m_edge, m_face))
        {
          m_face = null_face();
          // need to update face
          auto facesRange = incident_faces(m_edge);
          auto itF = facesRange.begin();

          for(; itF != facesRange.end(); ++itF)
          {
            edgesRange = (*itF)->GetIncidentEdges();
            auto itE = edgesRange.begin();

            for(; itE != edgesRange.end(); ++itE)
            {
              // we don't know in which order are the vertices in the AIF edge
              if(*itE == m_edge)
              {
                ++itE;
                if(itE == edgesRange.end())
                  itE = edgesRange.begin();

                if(common_vertex(m_edge, *itE) == get_target())
                  m_face = *itF;
                break;
              }
            }
          }
          if(m_face != null_face())
            edgesRange = m_face->GetIncidentEdges();
        }
        auto itE = edgesRange.begin();

        if(edgesRange.size() == 1)
        {
          if(*itE != m_edge)
          {
            found = true;
            edge = *itE;

            return AIFHalfEdge(get_target(),
                               opposite_vertex(edge, get_target()),
                               m_face,
                               edge,
                               false);
          }
        }
        else if(edgesRange.size() == 2)
        {
          if(*itE == m_edge)
            ++itE;

          found = true;
          edge = *itE;

          return AIFHalfEdge(get_target(),
                             opposite_vertex(edge, get_target()),
                             m_face,
                             edge,
                             false);
        }
        else
        {
          bool is_current_Edge_present = false;
          for(; itE != edgesRange.end(); ++itE)
          {
            if(*itE == m_edge)
            {
              is_current_Edge_present = true;
              continue;
            }
            auto verticesPair = (*itE)->GetIncidentVertices();
            vertex_type::ptr ptrV1 = verticesPair[0];
            vertex_type::ptr ptrV2 = verticesPair[1];

            // we don't know in which order are the vertices in the AIF edge
            if(ptrV1 == get_target() && ptrV2 != get_source())
            {
              if(found)
              {
                // edgeAmb = *itE;
                if(is_current_Edge_present)
                  found = false; // ambiguity detected
              }
              else
              {
                found = true;
                ptrVn = ptrV2;
                edge = *itE;
                // break;
              }
            }
            else if(ptrV2 == get_target() && ptrV1 != get_source())
            {
              if(found)
              {
                // edgeAmb = *itE;
                if(is_current_Edge_present)
                  found = false; // ambiguity detected
              }
              else
              {
                found = true;
                ptrVn = ptrV1;
                edge = *itE;
                // break;
              }
            }
          }

          if(!found)
          { // particular case due to CGAL low-level topological operations
            itE = edgesRange.begin();
            for(; itE != edgesRange.end(); ++itE)
            {
              auto verticesPair = (*itE)->GetIncidentVertices();
              vertex_type::ptr ptrV1 = verticesPair[0];
              vertex_type::ptr ptrV2 = verticesPair[1];

              // we don't know in which order are the vertices in the AIF edge
              if(*itE == m_edge)
              {
                ++itE; // we take the edge following the edge
                if(itE == edgesRange.end())
                  itE = edgesRange.begin();
                found = true;
                edge = *itE;
                ptrV1 = common_vertex(edge, m_edge);

                return AIFHalfEdge(
                    ptrV1, opposite_vertex(edge, ptrV1), m_face, edge, false);
              }
              else if(ptrV1 == ptrV2)
              {
                found = true;
                edge = *itE;
                return AIFHalfEdge(get_target(),
                                   opposite_vertex(edge, ptrV1),
                                   m_face,
                                   edge,
                                   false);
              }
            }
          }
        }
      }

      // abort if no edge with the right source is found
      if(!found)
        throw std::runtime_error(
            "Error in AIFHalfEdge::next(): next halfedge not found in face! "
            "Maybe the input halfedge is invalid (face/vertices mismatch).");

      // next edge found, return a halfedge
      return AIFHalfEdge(get_target(), ptrVn, m_face, edge, false);
    }

    /// prev of current halfedge
    AIFHalfEdge prev(const AIFMesh &m)
    {
      // parse edges around the face
      // look for the edge whose target is current halfedge source
      bool found = false;
      bool isBorderHalfedge = (m_face == null_face());
      edge_descriptor edge;
      vertex_type::ptr ptrVn;
      if(isBorderHalfedge)
      {
        // Check if target vertex is a 2-manifold vertex
        if(is_2_manifold_vertex(get_source()))
        { // in that case, we know there is only one border edge among
          // all edges incident to get_source()
          auto edgesRange = get_source()->GetIncidentEdges();
          auto itE = edgesRange.begin();
          for(; itE != edgesRange.end(); ++itE)
          {
            if(!is_surface_border_edge(*itE))
              continue;

            auto verticesPair = (*itE)->GetIncidentVertices();
            vertex_type::ptr ptrV1 = verticesPair[0];
            vertex_type::ptr ptrV2 = verticesPair[1];

            // we don't know in which order are the vertices in the AIF edge
            if(ptrV1 == get_source() && ptrV2 != get_target())
            {
              found = true;
              ptrVn = ptrV2;
              edge = *itE;
              break;
            }

            if(ptrV2 == get_source() && ptrV1 != get_target())
            {
              found = true;
              ptrVn = ptrV1;
              edge = *itE;
              break;
            }
          }
        }
        else
        {
          // We know that the previous edge belongs to another face if the
          // vertex is a cut-vertex
          if(is_cut_vertex(get_source()))
          {
            auto currentFaces = incident_faces(m_edge);
            auto facesRange = incident_faces(get_source());
            auto itF = facesRange.begin();
            int cpt = 0;
            face_descriptor rightFaceContainingPrevHalfEdge;
            for(; itF != facesRange.end(); ++itF)
            {
              if(*itF == currentFaces[0])
                continue;
              rightFaceContainingPrevHalfEdge = *itF;
              ++cpt;
            }
            if(cpt == 1)
            {
              auto edgesRange = incident_edges(rightFaceContainingPrevHalfEdge);
              auto itE = edgesRange.begin();
              edge_descriptor candidate1, candidate2;
              for(; itE != edgesRange.end(); ++itE)
              {
                if(((*itE)->get_first_vertex() == get_source()) ||
                   ((*itE)->get_second_vertex() == get_source()))
                {
                  if(candidate1 == null_edge())
                    candidate1 = *itE;
                  else if(candidate2 == null_edge())
                    candidate2 = *itE;
                  else
                    throw std::runtime_error(
                        "Error in AIFHalfEdge::prev(): prev halfedge cannot be "
                        "found in a non-manifold context. Get unkown case "
                        "while processing a cut-vertex.");
                }
              }
              found = true;
              if(is_surface_border_edge(candidate1))
              {
                if(is_surface_border_edge(candidate2))
                {
                  AIFHalfEdge h(candidate1->get_first_vertex(),
                                candidate1->get_second_vertex(),
                                rightFaceContainingPrevHalfEdge,
                                true);
                  if(h.prev(m).get_edge() == candidate2)
                  {
                    ptrVn = opposite_vertex(candidate2, get_source());
                    edge = candidate2;
                  }
                  else if(h.next(m).get_edge() == candidate2)
                  {
                    ptrVn = opposite_vertex(candidate1, get_source());
                    edge = candidate1;
                  }
                }
                else
                {
                  ptrVn = opposite_vertex(candidate1, get_source());
                  edge = candidate1;
                }
              }
              else if(is_surface_border_edge(candidate2))
              {
                ptrVn = opposite_vertex(candidate2, get_source());
                edge = candidate2;
              }
              else
              {
                found = false;
                throw std::runtime_error(
                    "Error in AIFHalfEdge::prev(): prev halfedge cannot be "
                    "found in a non-manifold context. Get unkown case while "
                    "processing a cut-vertex.");
              }
            }
            else
              throw std::runtime_error(
                  "Error in AIFHalfEdge::prev(): prev halfedge cannot be found "
                  "in a non-manifold context. Current vertex is a cut-vertex.");
          }

          if(is_degenerated_vertex(get_source()))
            throw std::runtime_error(
                "Error in AIFHalfEdge::prev(): prev halfedge cannot be found "
                "in a non-manifold context. Current vertex is degenerated.");

          if(is_isolated_vertex(get_source()))
            throw std::runtime_error(
                "Error in AIFHalfEdge::prev(): prev halfedge cannot be found "
                "in a non-manifold context. Current vertex is isolated.");

          throw std::runtime_error(
              "Error in AIFHalfEdge::prev(): prev halfedge cannot be found in "
              "a non-manifold context. An incident edge is non-manifold.");
        }
      }
      else
      {
        auto edgesRange = m_face->GetIncidentEdges();

        if((edgesRange.size() == 0) && !are_incident(m_edge, m_face))
        {
          m_face = null_face();
          // need to update face
          auto facesRange = incident_faces(m_edge);
          auto itF = facesRange.begin();

          for(; itF != facesRange.end(); ++itF)
          {
            edgesRange = (*itF)->GetIncidentEdges();
            auto itE = edgesRange.begin();

            for(; itE != edgesRange.end(); ++itE)
            {
              // we don't know in which order are the vertices in the AIF edge
              if(*itE == m_edge)
              {
                ++itE;
                if(itE == edgesRange.end())
                  itE = edgesRange.begin();

                if(common_vertex(m_edge, *itE) == get_target())
                  m_face = *itF;
                break;
              }
            }
          }
          if(m_face != null_face())
            edgesRange = m_face->GetIncidentEdges();
        }

        auto itE = edgesRange.begin();

        if(edgesRange.size() == 1)
        {
          if(*itE != m_edge)
          {
            found = true;
            edge = *itE;

            return AIFHalfEdge(
                opposite_vertex(edge, m_source), m_source, m_face, edge, false);
          }
        }
        else if(edgesRange.size() == 2)
        {
          if(*itE == m_edge)
            ++itE;

          found = true;
          edge = *itE;

          return AIFHalfEdge(
              opposite_vertex(edge, m_source), m_source, m_face, edge, false);
        }
        else
        {
          bool is_current_Edge_present = false;
          for(; itE != edgesRange.end(); ++itE)
          {
            if(*itE == m_edge)
            {
              is_current_Edge_present = true;
              continue;
            }
            auto verticesPair = (*itE)->GetIncidentVertices();
            vertex_type::ptr ptrV1 = verticesPair[0];
            vertex_type::ptr ptrV2 = verticesPair[1];

            // we don't know in which order are the vertices in the AIF edge
            if(ptrV1 == get_source() && ptrV2 != get_target())
            {
              if(found)
              {
                if(is_current_Edge_present)
                  found = false; // ambiguity detected
              }
              else
              {
                found = true;
                ptrVn = ptrV2;
                edge = *itE;
                // break;
              }
            }
            else if(ptrV2 == get_source() && ptrV1 != get_target())
            {
              if(found)
              {
                if(is_current_Edge_present)
                  found = false; // ambiguity detected
              }
              else
              {
                found = true;
                ptrVn = ptrV1;
                edge = *itE;
                // break;
              }
            }
          }
          if(!found)
          { // particular case due to CGAL low-level topological operations
            itE = edgesRange.begin();
            for(; itE != edgesRange.end(); ++itE)
            {
              auto verticesPair = (*itE)->GetIncidentVertices();
              vertex_type::ptr ptrV1 = verticesPair[0];
              vertex_type::ptr ptrV2 = verticesPair[1];

              // we don't know in which order are the vertices in the AIF edge
              if(*itE == m_edge)
              {
                if(itE == edgesRange.begin())
                  itE = edgesRange.end();
                --itE; // we take the edge preceding the edge

                found = true;
                edge = *itE;
                ptrV1 = common_vertex(edge, m_edge);
                return AIFHalfEdge(
                    opposite_vertex(edge, ptrV1), ptrV1, m_face, edge, false);
              }
              else if(ptrV1 == ptrV2)
              {
                found = true;
                edge = *itE;
                return AIFHalfEdge(get_target(),
                                   opposite_vertex(edge, ptrV1),
                                   m_face,
                                   edge,
                                   false);
              }
            }
          }
        }
      }

      // abort if no edge with the right source is found
      if(!found)
        throw std::runtime_error(
            "Error in AIFHalfEdge::prev(): prev halfedge not found in face! "
            "Maybe the input halfedge is invalid (face/vertices mismatch).");

      // prev edge found, return a halfedge
      return AIFHalfEdge(ptrVn, get_source(), m_face, edge, false);
    }

    /// opposite of current halfedge (may not be defined when edge is complex)
    AIFHalfEdge opposite(const AIFMesh &m)
    {
      edge_descriptor good_e = m_edge;
      face_descriptor op_f; // op_f can be null for border edge (not an issue)
      int cpt = 0;
      // parse edges around the face
      if(m_face == null_face())
      {
        auto facesRange = incident_faces(good_e);
        auto itF = facesRange.begin();

        for(; itF != facesRange.end(); ++itF, ++cpt)
        {
          op_f = *itF;
        }
        if(cpt != 1)
        {
          return AIFHalfEdge(
              get_target(),
              get_source(),
              null_face(),
              good_e,
              false); // for some Euler operations, we can handle non-valid
                      // halfedge at a time (e.g. split_vertex)
        }
        return AIFHalfEdge(get_target(), get_source(), op_f, good_e, false);
      }

      auto facesRange = incident_faces(good_e);
      auto itF = facesRange.begin();

      for(; itF != facesRange.end(); ++itF, ++cpt)
      {
        if(*itF != m_face)
          op_f = *itF;
      }
      if(cpt > 2)
      {
        throw std::invalid_argument(
            "AIFTopologyHelpers::opposite(h, m) -> several opposite edges "
            "exist, cannot proceed.");
      }

      return AIFHalfEdge(get_target(), get_source(), op_f, good_e, false);
    }

    /// source vertex getter
    inline vertex_descriptor get_source(void) const { return m_source; }
    /// target vertex getter
    inline vertex_descriptor get_target(void) const { return m_target; }
    /// face getter
    inline face_descriptor get_face(void) const { return m_face; }
    /// corresponding edge getter
    inline edge_descriptor get_edge(void) const { return m_edge; }

    /// source vertex setter
    void set_source(vertex_descriptor new_source) { m_source = new_source; }
    /// target vertex setter
    void set_target(vertex_descriptor new_target) { m_target = new_target; }
    /// face setter
    void set_face(face_descriptor new_face) { m_face = new_face; }
    /// edge setter
    void set_edge(edge_descriptor new_edge) { m_edge = new_edge; }

    /// ostream operator for debugging
    friend std::ostream &operator<<(std::ostream &stream, const AIFHalfEdge &h)
    {
      stream << "he(source=" << h.get_source() << ", target=" << h.get_target()
             << ", face=" << h.get_face() << ")";
      return stream;
    }

    // member attributes
  protected:
    vertex_descriptor m_source;
    vertex_descriptor
        m_target; // redundant information, but needed for some topological
                  // operations due to non persistency of halfedge

    face_descriptor m_face;
    edge_descriptor
        m_edge; // needed to handle modifications (the current
                // halfedge_descriptor is not a pointer or smart pointer, thus
                // all functions working on halfedge_descriptor works on copy
                // and not on original data). however the edge_descriptor is a
                // smart point, so we can work on original data.
  };

  typedef AIFHalfEdge halfedge_descriptor;

  /*!
   * 			Returns a null halfedge.
   * \return	The null halfedge.
   */
  static halfedge_descriptor null_halfedge() { return halfedge_descriptor(); }

  /*!
   * 			Return a halfedge incident to face.
   * \param	face	The face the halfedge must be incident to
   * \param	mesh	The involving mesh
   * \return	The halfedge incident to the face.
   */
  static halfedge_descriptor halfedge(face_descriptor face, const AIFMesh &mesh)
  {
    auto edgesRange = incident_edges(face);
    if(edgesRange.begin() == edgesRange.end())
      throw std::invalid_argument(
          "AIFTopologyHelpers::halfedge(f, m) -> no edge incident to face.");

    edge_descriptor e = *(edgesRange.begin());

    // the vertex that is common with the next edge in face-incident-edges
    // container must be the halfedge target
    edge_descriptor enext = *(edgesRange.begin() + 1);
    vertex_descriptor source, target;
    if(e->get_first_vertex() == enext->get_first_vertex() ||
       e->get_first_vertex() == enext->get_second_vertex())
    {
      target = e->get_first_vertex();
      source = e->get_second_vertex();
    }
    else if(e->get_second_vertex() == enext->get_first_vertex() ||
            e->get_second_vertex() == enext->get_second_vertex())
    {
      target = e->get_second_vertex();
      source = e->get_first_vertex();
    }
    else
      throw std::runtime_error(
          "AIFTopologyHelpers::halfedge(f, m) -> the second edge in "
          "face-incident-edges container don't share any vertex with the first "
          "edge in container.");

    // DBG*/ std::cout << __FILE__ << ":" << __LINE__ << " halfedge(face=" <<
    // face << ", mesh=" << &mesh << ")] = " << AIFHalfEdge(source, target, face)
    // << "    first_face_edge = " << e <<  "    second_face_edge = " << enext <<
    // std::endl;

    return AIFHalfEdge(source, target, face, e, false);
  }

  /*!
   * 			Return a halfedge incident to vertex.
   * \param	v	The vertex the halfedge must be incident to
   * \param	mesh	The involving mesh
   * \return	The halfedge incident to the face.
   */
  static halfedge_descriptor halfedge(vertex_descriptor v, const AIFMesh &mesh)
  {
    auto edgesRange = incident_edges(v);
    if(edgesRange.begin() == edgesRange.end())
      return null_halfedge(); // throw
                              // std::invalid_argument("AIFTopologyHelpers::halfedge(v,
                              // m) -> no edge incident to vertex.");

    auto itE = edgesRange.begin();
    edge_descriptor e = *itE;
    auto facesRange = incident_faces(e);
    while(itE != edgesRange.end())
    {
      if((facesRange.size() > 0) && (facesRange.size() < 3))
        break;
      ++itE;
      if(itE != edgesRange.end())
      {
        e = *itE;
        facesRange = incident_faces(e);
      }
    }
    if(itE == edgesRange.end())
      throw std::invalid_argument(
          "AIFTopologyHelpers::halfedge(v, m) -> cannot find a manifold edge.");

    if(facesRange.begin() == facesRange.end())
      throw std::invalid_argument("AIFTopologyHelpers::halfedge(v, m) -> no "
                                  "face incident to selected edge.");

    face_descriptor f = *(facesRange.begin());

    // select the right face
    auto edgesRangeF = incident_edges(f);
    auto iter = edgesRangeF.begin();
    edge_descriptor eF = *iter;
    ++iter;
    edge_descriptor eFnext = *iter;
    while(!(((eF->get_first_vertex() == v) &&
             (eF->get_second_vertex() == opposite_vertex(e, v))) ||
            ((eF->get_first_vertex() == opposite_vertex(e, v)) &&
             (eF->get_second_vertex() == v))))
    {
      eF = *iter;
      ++iter;
      if(iter == edgesRangeF.end())
        iter = edgesRangeF.begin();
      eFnext = *iter;
    }
    if(eF->get_first_vertex() == eFnext->get_first_vertex() ||
       eF->get_first_vertex() == eFnext->get_second_vertex())
    {
      if(eF->get_first_vertex() != v)
      {
        if(facesRange.size() == 1)
          f = null_face();
        else
          f = *(facesRange.begin() + 1);
      }
    }
    else if(eF->get_second_vertex() == eFnext->get_first_vertex() ||
            eF->get_second_vertex() == eFnext->get_second_vertex())
    {
      if(eF->get_second_vertex() != v)
      {
        if(facesRange.size() == 1)
          f = null_face();
        else
          f = *(facesRange.begin() + 1);
      }
    }
    ///////////////////////////////////////////////////////////////////
    return AIFHalfEdge(
        opposite_vertex(e, v),
        v,
        f,
        e,
        false); // here we cannot propose an halfedge outside the face!!
  }

  /*!
   * Return an halfedge src-->target if such an halfedge exist or
   * returns a null halfedge.
   *
   * \param  src     The source vertex of the * halfedge
   * \param  target  The target vertex of the halfedge (the halfedge is
   *                 pointing towards the target vertex)
   * \param  mesh    The involving mesh
   *
   * \return  The halfedge or a null halfedge.
   */
  static halfedge_descriptor
  halfedge(vertex_descriptor src, vertex_descriptor target, const AIFMesh &mesh)
  {
    edge_descriptor good_e = common_edge(src, target);
    if(good_e == null_edge())
      return halfedge_descriptor();
    ///////////////////////////////////////////////////////////////////
    face_descriptor op_f =
        null_face(); // op_f can be null for border edge (not an issue)
    ///////////////////////////////////////////////////////////////////
    int cpt = 0;
    auto facesRange = incident_faces(good_e);
    auto itF = facesRange.begin();

    for(; itF != facesRange.end(); ++itF, ++cpt)
    {
      edge_descriptor next_good_e;

      auto edgesRange = (*itF)->GetIncidentEdges();
      auto itE = edgesRange.begin();

      for(; itE != edgesRange.end(); ++itE)
      {
        if((((*itE)->get_first_vertex() == src) &&
            ((*itE)->get_second_vertex() == target)) ||
           (((*itE)->get_first_vertex() == target) &&
            ((*itE)->get_second_vertex() == src)))
        {
          good_e = *itE;
          ++itE;
          if(itE == edgesRange.end())
            itE = edgesRange.begin();
          next_good_e = *itE;
          // if the common vertex between good_e and next_good_e is src then,
          // the halfedge to return has no incident face. else, if the common
          // vertex between good_e and next_good_e is target, then the current
          // face is incident to halfedge.
          if(common_vertex(good_e, next_good_e) == target)
          {
            op_f = *itF;
          }
          break;
        }
      }
    }
    if(cpt == 0 || cpt > 2)
    {
      throw std::invalid_argument("AIFTopologyHelpers::halfedge(src, target, "
                                  "m) -> current halfedge is not manifold.");
    }
    ////////////////////////////////////////////////////////
    return halfedge_descriptor(src, target, op_f, good_e, false);
  }
  /*!
   * Return an halfedge corresponding to e edge if such an halfedge
   * exist or returns a null halfedge.
   *
   * \param  edge  The edge for which one halfedge is requested
   * \param  mesh  The involving mesh
   *
   * \return  The halfedge or a null halfedge.
   */
  static halfedge_descriptor halfedge(edge_descriptor edge, const AIFMesh &mesh)
  {
    if((edge->get_first_vertex() == null_vertex()) ||
       (edge->get_second_vertex() == null_vertex()))
    {
      halfedge_descriptor h = null_halfedge();
      h.set_source(edge->get_first_vertex());
      h.set_target(edge->get_second_vertex());
      h.set_edge(edge);
      return h;
    }
    else
      return halfedge(
          edge->get_first_vertex(), edge->get_second_vertex(), mesh);
  }
  /*!
   * Return an edge corresponding to h halfedge if such an halfedge
   * exist or returns a null edge.
   *
   * \param  h     The halfedge for which the corresponding edge is requested
   * \param  mesh  The involving mesh
   *
   * \return  The edge or a null edge.
   */
  static edge_descriptor edge(halfedge_descriptor h, const AIFMesh &mesh)
  {
    return h.get_edge();
  }
  ////////////////////////////////////////////////////////
  /*!
   * 			Set the target vertex of an halfedge. For MutableHalfedgeGraph
   * (CGAL).
   * \param	h	The halfedge for which the target vertex is to be set.
   * \param  target The vertex to use to set the halfedge target vertex.
   * \param	mesh	The involving mesh
   * \warning For the time being the new target modification is
   *          not taken into account because we work on a copy (not halfedge is
   * stored within the AIFMesh).
   */
  static void
  set_target(halfedge_descriptor h, vertex_descriptor target, AIFMesh &mesh)
  {
    edge_descriptor good_e = h.get_edge();
    if(good_e == null_edge())
      throw std::invalid_argument("AIFTopologyHelpers::set_target(h, v, m) -> "
                                  "halfedge does not belong to input mesh.");

    add_vertex_to_edge(
        good_e,
        target,
        vertex_position(
            good_e,
            opposite_vertex(good_e,
                            h.get_source()))); // h.get_source() is consistent,
                                               // not necesseraly h.get_target()
    h.set_target(
        target); // to update the halfedge (needed for following operations)
  }
  /*!
   * Set the incoming halfedge of vertex target to h. For
   * MutableHalfedgeGraph (CGAL).
   *
   * \param  target  The vertex to use to set the incoming halfedge.
   * \param  h       The halfedge.
   * \param  mesh    The involving mesh
   */
  static void
  set_halfedge(vertex_descriptor target, halfedge_descriptor h, AIFMesh &mesh)
  {
    edge_descriptor good_e = h.get_edge();
    if(good_e == null_edge())
      throw std::invalid_argument("AIFTopologyHelpers::set_halfedge(v, h, m) "
                                  "-> halfedge does not belong to input mesh.");

    if(std::find(target->GetIncidentEdges().begin(),
                 target->GetIncidentEdges().end(),
                 good_e) == target->GetIncidentEdges().end())
      add_edge_to_vertex(target, good_e);
  }
  /*!
   * Sets the successor of h1 around a face to h2, and the
   * prededecessor of h2 to h1. For MutableHalfedgeGraph (CGAL). Usage
   * precondition: h1 != h2 and h1 != opposite(h2, mesh)
   *
   * \param  h1    The halfedge for which the next halfedge is to be set.
   * \param  h2    The halfedge for which the previous halfedge is to be set.
   * \param  mesh  The involving mesh
   *
   * \warning  For the time being the h2.set_face(face) modification is not
   *           taken into account because we work on a copy (not halfedge
   *           is stored within the AIFMesh).
   */
  static void
  set_next(halfedge_descriptor h1, halfedge_descriptor h2, AIFMesh &mesh)
  {
    assert(h1 != h2);
    assert(h1 != h2.opposite(mesh));
    face_descriptor face = h1.get_face();
    if(face == null_face())
    {
      face = h2.get_face(); // some CGAL functions assume that a call to
                            // set_next can be done with h1 having a null face
                            // and h2 having a null face too
    }
    edge_descriptor prev_edge = h1.get_edge();
    if(prev_edge == null_edge())
      throw std::invalid_argument("AIFTopologyHelpers::set_next(h1, h2, m) -> "
                                  "halfedge h1 does not belong to input mesh.");

    edge_descriptor edge_to_insert = h2.get_edge();
    if(edge_to_insert == null_edge())
      throw std::invalid_argument("AIFTopologyHelpers::set_next(h1, h2, m) -> "
                                  "halfedge h2 does not belong to input mesh.");
    halfedge_descriptor current_next_h1 = h1.next(mesh),
                        current_next_h1_op = current_next_h1.opposite(mesh);
    if(current_next_h1 != h2)
    {
      if(face != null_face())
        remove_edge_from_face(face, current_next_h1.get_edge());
      if(are_incident(common_vertex(h1.get_edge(), current_next_h1.get_edge()),
                      current_next_h1.get_edge()))
        remove_edge_from_vertex(
            common_vertex(h1.get_edge(), current_next_h1.get_edge()),
            current_next_h1.get_edge());
      if(face != h2.get_face())
      {
        std::vector< edge_descriptor > edges_to_insert_in_face;
        if(face->GetDegree() > 0)
        {
          halfedge_descriptor tmp = h2;
          while(tmp != current_next_h1_op)
          {
            edges_to_insert_in_face.push_back(tmp.get_edge());
            tmp = tmp.next(mesh);
          }

          if(!are_incident(edge_to_insert, face))
          {
            link_edges_and_face_around_face_after_edge(
                edges_to_insert_in_face, face, prev_edge);
            unlink_edges_and_face(
                edges_to_insert_in_face,
                current_next_h1_op.get_face()); // to ensure 2-manifoldness
          }
          h2.set_face(
              face); // to update the halfedge (needed for following operations)
        }
        else
        {
          if(is_dangling_edge(current_next_h1.get_edge()))
          {
            if(std::find(current_next_h1.get_edge()
                             ->get_first_vertex()
                             ->GetIncidentEdges()
                             .begin(),
                         current_next_h1.get_edge()
                             ->get_first_vertex()
                             ->GetIncidentEdges()
                             .end(),
                         current_next_h1.get_edge()) !=
               current_next_h1.get_edge()
                   ->get_first_vertex()
                   ->GetIncidentEdges()
                   .end())
              remove_edge_from_vertex(
                  current_next_h1.get_edge()->get_first_vertex(),
                  current_next_h1.get_edge());

            if(std::find(current_next_h1.get_edge()
                             ->get_second_vertex()
                             ->GetIncidentEdges()
                             .begin(),
                         current_next_h1.get_edge()
                             ->get_second_vertex()
                             ->GetIncidentEdges()
                             .end(),
                         current_next_h1.get_edge()) !=
               current_next_h1.get_edge()
                   ->get_second_vertex()
                   ->GetIncidentEdges()
                   .end())
              remove_edge_from_vertex(
                  current_next_h1.get_edge()->get_second_vertex(),
                  current_next_h1.get_edge());
          }
          face_descriptor common_face;
          auto h1FaceRange = incident_faces(h1.get_edge());
          auto h2FaceRange = incident_faces(h2.get_edge());
          auto h1itF = h1FaceRange.begin();
          while(h1itF != h1FaceRange.end())
          {
            auto h2itF = h2FaceRange.begin();
            while(h2itF != h2FaceRange.end())
            {
              if(*h1itF == *h2itF)
                common_face = *h1itF;
              ++h2itF;
            }
            ++h1itF;
          }

          h1.set_face(common_face);
          h2.set_face(common_face);
        }
      }
      // else
      //  throw std::runtime_error("AIFTopologyHelpers::set_next(h1, h2, m) ->
      //  case next(h1,m)!=h2 and face(h1,m)==face(h2,m) not yet managed.");
    }
    assert(h1.next(mesh) == h2);
    assert(h2.prev(mesh) == h1);
  }

  /*!
   * 			Sets the corresponding face of h to face if h and face were not
   * incident. For MutableHalfedgeGraph (CGAL). Usage precondition: edge(h,
   * mesh) != null_edge()
   * \param	h	The halfedge for which the face is to be set.
   * \param	face	The face to use.
   * \param	mesh	The involving mesh
   * \warning For the time being the h2.set_face(face) modification is
   *          not taken into account because we work on a copy (not halfedge is
   * stored within the AIFMesh).
   */
  static void
  set_face(halfedge_descriptor h, face_descriptor face, AIFMesh &mesh)
  {
    edge_descriptor edge = h.get_edge();
    if(edge == null_edge())
      throw std::invalid_argument("AIFTopologyHelpers::set_face(h, f, m) -> "
                                  "halfedge does not belong to input mesh.");

    if((face != null_face()) && !are_incident(edge, face))
    {
      add_face_to_edge(edge, face);
    }
    else if((face == null_face()) && (h.get_face() != null_face()))
    {
      if(are_incident(edge, h.get_face()))
        unlink_edge_and_face(edge, h.get_face());
    }
    if(h.get_face() != face)
      h.set_face(
          face); // to update the halfedge (needed for following operations)
  }
  /*!
   * Sets the corresponding halfedge of face to h if face and h were
   * not incident. For MutableHalfedgeGraph (CGAL). Usage precondition: face !=
   * null_face() and edge(h, mesh) != null_edge()
   *
   * \param  h     The halfedge for which the face is to be set.
   * \param  face  The face to use.
   * \param  mesh  The involving mesh
   *
   * \warning  For the time being the h2.set_face(face) modification is
   *           not taken into account because we work on a copy (not halfedge
   *           is stored within the AIFMesh).
   */
  static void
  set_halfedge(face_descriptor face, halfedge_descriptor h, AIFMesh &mesh)
  {
    if(face == null_face())
      throw std::invalid_argument(
          "AIFTopologyHelpers::set_halfedge(f, h, m) -> face is a null face.");
    edge_descriptor edge = h.get_edge();
    if(edge == null_edge())
      throw std::invalid_argument("AIFTopologyHelpers::set_halfedge(f, h, m) "
                                  "-> halfedge does not belong to input mesh.");
    if(!are_incident(face, edge))
    {
      add_edge_to_face(face, edge);
    }
  }
  // END OF HALFEDGE STUFF
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
};

} // namespace AIF
} // namespace DataStructures
} // namespace FEVV
