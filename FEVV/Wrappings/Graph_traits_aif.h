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

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <boost/iterator/transform_iterator.hpp>
//#include "FEVV/Wrappings/Graph_properties_aif.h" // Do not include it!

#include <algorithm>

#if defined(BOOST_MSVC)
#pragma warning(push)
#pragma warning(disable : 4267)
#endif


// shortcut for mesh type
typedef FEVV::DataStructures::AIF::AIFMesh AIFMeshT;
typedef FEVV::DataStructures::AIF::AIFTopologyHelpers AIFHelpers;


namespace boost {
template<>
struct vertex_property_type< FEVV::DataStructures::AIF::AIFMesh >
{
  // typedef AIFMeshT::vertex_type::ptr type;
  typedef AIFMeshT::Point type;
};

template<>
struct vertex_property_type< const FEVV::DataStructures::AIF::AIFMesh >
{
  // typedef AIFMeshT::vertex_type::ptr type;
  typedef AIFMeshT::Point type;
};

template<>
struct edge_property_type< FEVV::DataStructures::AIF::AIFMesh >
{
  typedef void type;
};

template<>
struct edge_property_type< const FEVV::DataStructures::AIF::AIFMesh >
{
  typedef void type;
};

template<>
struct graph_traits< FEVV::DataStructures::AIF::AIFMesh >
{
public:
  // Graph
  typedef AIFMeshT::vertex_type::ptr vertex_descriptor;
  typedef AIFMeshT::Point vertex_property_type;
  typedef AIFMeshT::edge_type::ptr edge_descriptor;
  typedef void edge_property_type;
  typedef boost::undirected_tag directed_category;
  typedef boost::disallow_parallel_edge_tag edge_parallel_category;
  // typedef boost::bidirectional_graph_tag		   traversal_category;
  struct my_traversal_category : public boost::vertex_list_graph_tag,
                                 public boost::edge_list_graph_tag,
                                 public boost::adjacency_graph_tag,
                                 // public boost::incidence_graph_tag,
                                 public boost::bidirectional_graph_tag
  {
  };
  typedef my_traversal_category traversal_category;

  // FaceGraph
  typedef AIFMeshT::face_type::ptr face_descriptor;

  // FaceListGraph
  typedef AIFMeshT::FaceContainerType::iterator face_iterator;
  typedef AIFHelpers::size_type faces_size_type;

  // HalfedgeGraph
  typedef AIFHelpers::AIFHalfEdge halfedge_descriptor;
  // typedef AIFMeshT::EdgeContainerType::iterator	   halfedge_iterator;

  // VertexListGraph
  typedef AIFMeshT::VertexContainerType::iterator vertex_iterator;
  typedef AIFHelpers::size_type vertices_size_type;

  // EdgeListGraph
  typedef AIFMeshT::EdgeContainerType::iterator edge_iterator;
  typedef AIFHelpers::size_type edges_size_type;

  // BGL valid expression:
  // http://www.boost.org/doc/libs/1_55_0/libs/graph/doc/Graph.html
  static inline vertex_descriptor null_vertex() { return vertex_descriptor(); }

  // CGAL-BGL HalfedgeGraph Concept
  static inline halfedge_descriptor null_halfedge()
  {
    return halfedge_descriptor();
  }

  // No Concept
  static inline edge_descriptor null_edge() { return edge_descriptor(); }

  // CGAL-BGL FaceGraph Concept
  static inline face_descriptor null_face() { return face_descriptor(); }

  // IncidenceGraph Concept
  typedef AIFHelpers::degree_size_type degree_size_type;
  typedef AIFHelpers::adjacency_iterator adjacency_iterator;

  template< typename CellType = vertex_descriptor >
  struct IncidenceTraits // default is the 1D case
  {
    typedef AIFHelpers::out_edge_iterator out_edge_iterator;
    typedef AIFHelpers::out_edge_iterator in_edge_iterator;
  };

  typedef IncidenceTraits< vertex_descriptor >::out_edge_iterator
      out_edge_iterator;
  typedef IncidenceTraits< vertex_descriptor >::in_edge_iterator
      in_edge_iterator;
};

template<>
struct graph_traits< FEVV::DataStructures::AIF::AIFMesh >::IncidenceTraits<
    typename boost::graph_traits<
        FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor >
{
  // MUTS FaceIncidentGraph Concept
  typedef AIFHelpers::incident_face_iterator in_edge_iterator;
};
template<>
struct graph_traits< FEVV::DataStructures::AIF::AIFMesh >::IncidenceTraits<
    typename boost::graph_traits<
        FEVV::DataStructures::AIF::AIFMesh >::face_descriptor >
{
  // MUTS FaceIncidentGraph Concept
  typedef AIFHelpers::incident_edge_iterator out_edge_iterator;
};

template<>
struct graph_traits< const FEVV::DataStructures::AIF::AIFMesh >
    : public graph_traits< FEVV::DataStructures::AIF::AIFMesh >
{
};

} // namespace boost

namespace FEVV {
namespace DataStructures {
namespace AIF {

// Essential free functions specialization for AIF.

// See http://www.boost.org/doc/libs/1_61_0/libs/graph/doc/graph_concepts.html
// for BGL concepts description.

// See http://doc.cgal.org/latest/BGL/group__PkgBGLConcepts.html
// for CGAL-BGL concepts description.

// BGL VertexListGraph
//!
//! \brief  Returns the iterator range of the vertices of the mesh.
//!
inline std::pair< typename boost::graph_traits<
                      FEVV::DataStructures::AIF::AIFMesh >::vertex_iterator,
                  typename boost::graph_traits<
                      FEVV::DataStructures::AIF::AIFMesh >::vertex_iterator >
vertices(const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  // TODO-elo Fix const_cast when CGAL will be const correct
  auto vertices_range =
      (const_cast< FEVV::DataStructures::AIF::AIFMesh & >(sm)).GetVertices();

  return std::pair< typename boost::graph_traits<
                        FEVV::DataStructures::AIF::AIFMesh >::vertex_iterator,
                    typename boost::graph_traits<
                        FEVV::DataStructures::AIF::AIFMesh >::vertex_iterator >(
      vertices_range.begin(), vertices_range.end());
}

// BGL VertexListGraph
//!
//! \brief  Returns an upper bound of the number of vertices of the mesh.
//!
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::vertices_size_type
num_vertices(const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  return static_cast< typename boost::graph_traits<
      FEVV::DataStructures::AIF::AIFMesh >::vertices_size_type >(
      sm.GetNumberOfVertices());
}

// BGL VertexIndexGraph Concept
inline void
renumber_vertex_indices(const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  auto vertex_range_pair = vertices(sm);
  auto vertex_iter = vertex_range_pair.first;
  int index = 0;
  for(; vertex_iter != vertex_range_pair.second; ++vertex_iter)
  {
    (*vertex_iter)->SetIndex(index++);
  }
}


// BGL AdjacencyGraph
inline std::pair< typename boost::graph_traits<
                      FEVV::DataStructures::AIF::AIFMesh >::adjacency_iterator,
                  typename boost::graph_traits<
                      FEVV::DataStructures::AIF::AIFMesh >::adjacency_iterator >
adjacent_vertices(typename boost::graph_traits<
                      FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor v,
                  const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  auto vertices_range = AIFHelpers::adjacent_vertices(v);
  return std::pair<
      typename boost::graph_traits<
          FEVV::DataStructures::AIF::AIFMesh >::adjacency_iterator,
      typename boost::graph_traits<
          FEVV::DataStructures::AIF::AIFMesh >::adjacency_iterator >(
      vertices_range.begin(), vertices_range.end());
}

// BGL EdgeListGraph
//!
//! \brief  Returns the iterator range of the edges of the mesh.
//!
inline std::pair< typename boost::graph_traits<
                      FEVV::DataStructures::AIF::AIFMesh >::edge_iterator,
                  typename boost::graph_traits<
                      FEVV::DataStructures::AIF::AIFMesh >::edge_iterator >
edges(const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  // TODO-elo Fix const_cast when CGAL will be const correct
  auto edges_range =
      (const_cast< FEVV::DataStructures::AIF::AIFMesh & >(sm)).GetEdges();

  return std::pair< typename boost::graph_traits<
                        FEVV::DataStructures::AIF::AIFMesh >::edge_iterator,
                    typename boost::graph_traits<
                        FEVV::DataStructures::AIF::AIFMesh >::edge_iterator >(
      edges_range.begin(), edges_range.end());
}

// BGL EdgeListGraph
//!
//! \brief  Returns an upper bound of the number of edges of the graph.
//!
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::edges_size_type
num_edges(const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  return static_cast< typename boost::graph_traits<
      FEVV::DataStructures::AIF::AIFMesh >::edges_size_type >(
      sm.GetNumberOfEdges());
}

// BGL EdgeIndexGraph Concept 
inline
void
renumber_edge_indices(const FEVV::DataStructures::AIF::AIFMesh& sm)
{
  auto edgeRangePair = edges(sm);
  auto edgeIter = edgeRangePair.first;
  int index = 0;
  for (; edgeIter != edgeRangePair.second; ++edgeIter) {
    (*edgeIter)->SetIndex(index++);
  }
}

// CGAL FaceListGraph
//!
//! \brief  Returns an upper bound of the number of faces of the graph.
//!
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::faces_size_type
num_faces(const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  return static_cast< typename boost::graph_traits<
      FEVV::DataStructures::AIF::AIFMesh >::faces_size_type >(
      sm.GetNumberOfFaces());
}

// CGAL-BGL HalfedgeGraph Concept
//!
//! \brief  Returns a halfedge with target v.
//!
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::halfedge_descriptor
halfedge(typename boost::graph_traits<
             FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor v,
         const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  return AIFHelpers::halfedge(v, sm);
}

// CGAL-BGL HalfedgeGraph Concept
//!
//! \brief  Returns the halfedge with source u and target v. The Boolean is true
//!         if this halfedge exists.
//!
inline std::pair< typename boost::graph_traits<
                      FEVV::DataStructures::AIF::AIFMesh >::halfedge_descriptor,
                  bool >
halfedge(
    typename boost::graph_traits<
        FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor u, /// source
    typename boost::graph_traits<
        FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor v, /// target
    const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  typedef typename boost::graph_traits< FEVV::DataStructures::AIF::AIFMesh >
      GraphTrait;
  typename GraphTrait::halfedge_descriptor h = AIFHelpers::halfedge(u, v, sm);

  return std::pair< typename GraphTrait::halfedge_descriptor, bool >(
      h, (h != GraphTrait::null_halfedge()));
}

// CGAL-BGL HalfedgeGraph Concept
//!
//! \brief  Returns one of the halfedges corresponding to e.
//!
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::halfedge_descriptor
halfedge(typename boost::graph_traits<
             FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor e,
         const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  return AIFHelpers::halfedge(e, sm);
}

// CGAL-BGL HalfedgeGraph Concept
//!
//! \brief  Returns the edge corresponding to h and opposite(h).
//!
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor
edge(typename boost::graph_traits<
         FEVV::DataStructures::AIF::AIFMesh >::halfedge_descriptor h,
     const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  return AIFHelpers::edge(h, sm);
}

// BGL AdjacencyMatrix Concept
// The edge() function return in constant time.
//!
//! \brief  Returns the edge with extremities u and v.
//!
inline std::pair< typename boost::graph_traits<
                      FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor,
                  bool >
edge(typename boost::graph_traits<
         FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor u,
     typename boost::graph_traits<
         FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor v,
     const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  typedef typename boost::graph_traits< FEVV::DataStructures::AIF::AIFMesh >
      GraphTrait;
  typename GraphTrait::edge_descriptor e = AIFHelpers::common_edge(u, v);

  if(e != GraphTrait::null_edge())
  {
    if((e->get_first_vertex() != u) &&
       (e->get_second_vertex() != v)) // to get the rigth edge on demand (just
                                      // swapping its vertices when needed)
      AIFHelpers::swap_vertices(e);
  }
  return std::pair< typename GraphTrait::edge_descriptor, bool >(
      e, (e != GraphTrait::null_edge()));
}

// BGL IncidenceGraph Concept
//!
//! \brief  Returns the source vertex of e.
//!
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor
source(typename boost::graph_traits<
           FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor e,
       const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  return e->get_first_vertex();
}

// BGL IncidenceGraph Concept
//!
//! \brief  Returns the target vertex of e.
//!
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor
target(typename boost::graph_traits<
           FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor e,
       const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  return e->get_second_vertex();
}

// BGL IncidenceGraph Concept
///  The source vertex of an edge obtained via an out edge iterator is
///  guaranteed (for both directed and undirected graphs) to be the vertex u
///  used in the call to out_edges(u, g) and the target vertex must be a vertex
///  adjacent to u
inline std::pair< typename boost::graph_traits<
                      FEVV::DataStructures::AIF::AIFMesh >::out_edge_iterator,
                  typename boost::graph_traits<
                      FEVV::DataStructures::AIF::AIFMesh >::out_edge_iterator >
out_edges(typename boost::graph_traits<
              FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor u,
          const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  auto edges_range = u->GetIncidentEdges();
  auto e_it = edges_range.begin(), e_ite = edges_range.end();
  // Issue: this method is not constant as boost concept requires
  // (http://www.boost.org/doc/libs/1_61_0/libs/graph/doc/IncidenceGraph.html),
  // it is linear in the number of out-edges.
  // But making it constant would imply commenting the next 3 lines of source
  // codes, which would imply that each edge added to the vector of incidence
  // edges of u is oriented outward u. Thus for each edge (u, v) u and v would
  // need to store 2 edges in opposite direction, but we do not want to use that
  // extra memory.
  for(; e_it != e_ite; ++e_it)
    if((*e_it)->get_first_vertex() !=
       u) // get_first_vertex() is supposed to be the source vertex
      AIFHelpers::swap_vertices(*e_it);

  return std::pair<
      typename boost::graph_traits<
          FEVV::DataStructures::AIF::AIFMesh >::out_edge_iterator,
      typename boost::graph_traits<
          FEVV::DataStructures::AIF::AIFMesh >::out_edge_iterator >(
      edges_range.begin(), edges_range.end());
}

// BGL BidirectionalGraph Concept
///  For both directed and undirected graphs, the target of an out-edge is
///  required to be vertex u and the source is required to be a vertex that is
///  adjacent to u.
inline std::pair< typename boost::graph_traits<
                      FEVV::DataStructures::AIF::AIFMesh >::in_edge_iterator,
                  typename boost::graph_traits<
                      FEVV::DataStructures::AIF::AIFMesh >::in_edge_iterator >
in_edges(typename boost::graph_traits<
             FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor u,
         const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  auto edges_range = u->GetIncidentEdges();
  auto e_it = edges_range.begin(), e_ite = edges_range.end();
  for(; e_it != e_ite; ++e_it)
    if((*e_it)->get_second_vertex() !=
       u) // get_second_vertex() is supposed to be the target vertex
      AIFHelpers::swap_vertices(*e_it);

  return std::pair<
      typename boost::graph_traits<
          FEVV::DataStructures::AIF::AIFMesh >::in_edge_iterator,
      typename boost::graph_traits<
          FEVV::DataStructures::AIF::AIFMesh >::in_edge_iterator >(
      edges_range.begin(), edges_range.end());
}

// MUTS FaceIncidentGraph concept refined by CellIncidenceGraph Concept [see
// out_edges in CellIncidenceGraphConceptPage]
inline std::pair<
    typename boost::graph_traits< FEVV::DataStructures::AIF::AIFMesh >::
        IncidenceTraits< typename boost::graph_traits<
            FEVV::DataStructures::AIF::AIFMesh >::face_descriptor >::
            out_edge_iterator,
    typename boost::graph_traits< FEVV::DataStructures::AIF::AIFMesh >::
        IncidenceTraits< typename boost::graph_traits<
            FEVV::DataStructures::AIF::AIFMesh >::face_descriptor >::
            out_edge_iterator >
out_edges(typename boost::graph_traits<
              FEVV::DataStructures::AIF::AIFMesh >::face_descriptor f,
          const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  auto edge_range = AIFHelpers::incident_edges(f);

  return std::pair<
      typename boost::graph_traits< FEVV::DataStructures::AIF::AIFMesh >::
          IncidenceTraits< typename boost::graph_traits<
              FEVV::DataStructures::AIF::AIFMesh >::face_descriptor >::
              out_edge_iterator,
      typename boost::graph_traits< FEVV::DataStructures::AIF::AIFMesh >::
          IncidenceTraits< typename boost::graph_traits<
              FEVV::DataStructures::AIF::AIFMesh >::face_descriptor >::
              out_edge_iterator >(edge_range.begin(), edge_range.end());
}
// MUTS FaceIncidentGraph concept refined by CellIncidenceGraph Concept [see
// in_edges in CellIncidenceGraphConceptPage]
inline std::pair<
    typename boost::graph_traits< FEVV::DataStructures::AIF::AIFMesh >::
        IncidenceTraits< typename boost::graph_traits<
            FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor >::
            in_edge_iterator,
    typename boost::graph_traits< FEVV::DataStructures::AIF::AIFMesh >::
        IncidenceTraits< typename boost::graph_traits<
            FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor >::
            in_edge_iterator >
in_edges(typename boost::graph_traits<
             FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor e,
         const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  auto face_range = AIFHelpers::incident_faces(e);

  return std::pair<
      typename boost::graph_traits< FEVV::DataStructures::AIF::AIFMesh >::
          IncidenceTraits< typename boost::graph_traits<
              FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor >::
              in_edge_iterator,
      typename boost::graph_traits< FEVV::DataStructures::AIF::AIFMesh >::
          IncidenceTraits< typename boost::graph_traits<
              FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor >::
              in_edge_iterator >(face_range.begin(), face_range.end());
}

// BGL IncidenceGraph Concept
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::degree_size_type
out_degree(typename boost::graph_traits<
               FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor u,
           const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  return u->GetIncidentEdges().size();
}

// BGL BidirectionalGraph Concept
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::degree_size_type
in_degree(typename boost::graph_traits<
              FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor u,
          const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  return u->GetIncidentEdges().size();
}

// BGL BidirectionalGraph Concept
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::degree_size_type
degree(typename boost::graph_traits<
           FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor u,
       const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  return u->GetIncidentEdges().size();
}

// MUTS FaceIncidentGraph concept refined by CellIncidenceGraph Concept
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::degree_size_type
degree(typename boost::graph_traits<
           FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor e,
       const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  auto face_range = AIFHelpers::incident_faces(e);

  return face_range.size();
}

// CGALFaceGraph Concept
//!
//! \brief  Returns the number of halfedges incident to face f.
//!
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::degree_size_type
degree(typename boost::graph_traits<
           FEVV::DataStructures::AIF::AIFMesh >::face_descriptor f,
       const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  return f->GetDegree();
}

// CGAL-BGL HalfedgeGraph Concept
//!
//! \brief  Returns the source vertex of h.
//!
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor
source(typename boost::graph_traits<
           FEVV::DataStructures::AIF::AIFMesh >::halfedge_descriptor h,
       const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  return h.get_source();
}

// CGAL-BGL HalfedgeGraph Concept
//!
//! \brief  Returns the next halfedge around its face.
//!
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::halfedge_descriptor
next(typename boost::graph_traits<
         FEVV::DataStructures::AIF::AIFMesh >::halfedge_descriptor h,
     const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  return h.next(sm);
}

// CGAL-BGL HalfedgeGraph Concept
//!
//! \brief  Returns the previous halfedge around its face.
//!
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::halfedge_descriptor
prev(typename boost::graph_traits<
         FEVV::DataStructures::AIF::AIFMesh >::halfedge_descriptor h,
     const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  return h.prev(sm);
}

// CGAL-BGL HalfedgeGraph Concept
//!
//! \brief  Returns the halfedge with source and target swapped.
//!
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::halfedge_descriptor
opposite(typename boost::graph_traits<
             FEVV::DataStructures::AIF::AIFMesh >::halfedge_descriptor h,
         const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  return h.opposite(sm);
}

// CGAL-BGL HalfedgeGraph Concept
//!
//! \brief  Returns the target vertex of h.
//!
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor
target(typename boost::graph_traits<
           FEVV::DataStructures::AIF::AIFMesh >::halfedge_descriptor h,
       const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  return h.get_target();
}

// CGAL-BGL FaceGraph Concept
//!
//! \brief  Returns a halfedge incident to face f.
//!
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::halfedge_descriptor
halfedge(typename boost::graph_traits<
             FEVV::DataStructures::AIF::AIFMesh >::face_descriptor f,
         const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  return AIFHelpers::halfedge(f, sm);
}

// CGAL-BGL FaceGraph Concept
//!
//! \brief  Returns the face incident to halfedge h.
//!
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::face_descriptor
face(typename boost::graph_traits<
         FEVV::DataStructures::AIF::AIFMesh >::halfedge_descriptor h,
     const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  return h.get_face();
}

// CGAL-BGL FaceListGraph Concept
//!
//! \brief  Returns an iterator range over all faces of the mesh.
//!
inline std::pair< typename boost::graph_traits<
                      FEVV::DataStructures::AIF::AIFMesh >::face_iterator,
                  typename boost::graph_traits<
                      FEVV::DataStructures::AIF::AIFMesh >::face_iterator >
faces(const FEVV::DataStructures::AIF::AIFMesh &sm)
{
  // TODO-elo Fix const_cast when CGAL will be const correct
  auto faces_range =
      (const_cast< FEVV::DataStructures::AIF::AIFMesh & >(sm)).GetFaces();

  return std::pair< typename boost::graph_traits<
                        FEVV::DataStructures::AIF::AIFMesh >::face_iterator,
                    typename boost::graph_traits<
                        FEVV::DataStructures::AIF::AIFMesh >::face_iterator >(
      faces_range.begin(), faces_range.end());
}

// BGL VertexMutableGraph Concept
//!
//! \brief  Adds a new vertex to the graph without initializing the
//! connectivity.
//!
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor
add_vertex(FEVV::DataStructures::AIF::AIFMesh &sm)
{
  return AIFHelpers::add_vertex(sm);
}

// CGAL MutableHalfedgeGraph Concept
inline void
add_vertex(typename boost::graph_traits<
               FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor v,
           FEVV::DataStructures::AIF::AIFMesh &sm)
{
  typedef typename boost::graph_traits< FEVV::DataStructures::AIF::AIFMesh >
      GraphTrait;
  if(v == GraphTrait::null_vertex())
  {
    throw std::invalid_argument(
        "Helpers::add_vertex(v, m) -> vertex is null, cannot add it to mesh.");
  }
  else
  {
    AIFHelpers::add_vertex(v, sm);
  }
}

// CGAL MutableHalfedgeGraph Concept
//!
//! \brief  Removes v from the mesh.
//!
inline void
remove_vertex(typename boost::graph_traits<
                  FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor v,
              FEVV::DataStructures::AIF::AIFMesh &sm)
{
  AIFHelpers::remove_vertex(v, sm);
}

// CGAL MutableHalfedgeGraph Concept
//!
//! \brief  Sets the target vertex of h and the source of opposite(h) to v.
//!
inline void
set_target(typename boost::graph_traits< FEVV::DataStructures::AIF::AIFMesh >::
               halfedge_descriptor /*&*/ h, // cannot use a ref
           typename boost::graph_traits<
               FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor v,
           FEVV::DataStructures::AIF::AIFMesh &sm)
{
  AIFHelpers::set_target(h, v, sm);
}
/*
inline
void set_source(typename
boost::graph_traits<FEVV::DataStructures::AIF::AIFMesh>::edge_descriptor e,
  typename
boost::graph_traits<FEVV::DataStructures::AIF::AIFMesh>::vertex_descriptor v,
  FEVV::DataStructures::AIF::AIFMesh& sm)
{
  AIFHelpers::link_vertex_and_edge(v, e, AIFHelpers::vertex_position(e,
source(e, sm)));
}
inline
void set_target(typename
boost::graph_traits<FEVV::DataStructures::AIF::AIFMesh>::edge_descriptor e,
                typename
boost::graph_traits<FEVV::DataStructures::AIF::AIFMesh>::vertex_descriptor v,
                FEVV::DataStructures::AIF::AIFMesh& sm)
{
  AIFHelpers::link_vertex_and_edge(v, e, AIFHelpers::vertex_position(e,
target(e, sm)));
}*/

// CGAL MutableHalfedgeGraph Concept
//!
//! \brief  Sets the halfedge of v to h. The target vertex of h must be v.
//!
inline void
set_halfedge(typename boost::graph_traits<
                 FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor v,
             typename boost::graph_traits<
                 FEVV::DataStructures::AIF::AIFMesh >::halfedge_descriptor h,
             FEVV::DataStructures::AIF::AIFMesh &sm)
{
  AIFHelpers::set_halfedge(v, h, sm);
}

// CGAL MutableHalfedgeGraph Concept
//!
//! \brief  Sets the successor of h1 around a face to h2, and the prededecessor
//! of h2 to h1.
//!
inline void
set_next(typename boost::graph_traits<
             FEVV::DataStructures::AIF::AIFMesh >::halfedge_descriptor h1,
         typename boost::graph_traits<
             FEVV::DataStructures::AIF::AIFMesh >::halfedge_descriptor h2,
         FEVV::DataStructures::AIF::AIFMesh &sm)
{
  AIFHelpers::set_next(h1, h2, sm);
}
//#define USE_ADD_EDGE_AIF 
#ifdef USE_ADD_EDGE_AIF // conflict with function in Euler_operations.h (CGAL)
/// BGL MutableEdgeGraph Concept
/// Inserts the edge (u,v) into the graph, and returns an edge descriptor
/// pointing to the new edge. If the graph disallows parallel edges, and the
/// edge (u,v) is already in the graph, then the bool flag returned is false and
/// the returned edge descriptor points to the already existing edge. Possible
/// conflict with a function with tha same name and signature in
/// Euler_operations.h.
inline std::pair< typename boost::graph_traits<
                      FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor,
                  bool >
add_edge(typename boost::graph_traits<
             FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor u,
         typename boost::graph_traits<
             FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor v,
         FEVV::DataStructures::AIF::AIFMesh &sm)
{
  return AIFHelpers::add_edge(u, v, sm);
}
#endif

// CGAL MutableHalfedgeGraph Concept
//!
//! \brief  Adds two opposite halfedges to the graph without initializing the
//! connectivity.
//!
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor
add_edge(FEVV::DataStructures::AIF::AIFMesh &sm)
{
  return AIFHelpers::add_edge(sm);
}


// NOT IN MutableHalfedgeGraph Concept
inline void
add_edge(typename boost::graph_traits<
             FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor e,
         FEVV::DataStructures::AIF::AIFMesh &sm)
{
  typedef typename boost::graph_traits< FEVV::DataStructures::AIF::AIFMesh >
      GraphTrait;
  if(e == GraphTrait::null_edge())
  {
    throw std::invalid_argument(
        "Helpers::add_edge(e, m) -> edge is null, cannot add it to mesh.");
  }
  else
  {
    AIFHelpers::add_edge(e, sm);
  }
}


// CGAL MutableHalfedgeGraph Concept
//!
//! \brief  Remove edge e.
//!
inline void
remove_edge(typename boost::graph_traits<
                FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor e,
            FEVV::DataStructures::AIF::AIFMesh &sm)
{
  AIFHelpers::remove_edge(e, sm);
}

// BGL MutableEdgeGraph Concept
//!
//! \brief  Remove the edge with extremities u and v.
//!
inline void
remove_edge(typename boost::graph_traits<
                FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor u,
            typename boost::graph_traits<
                FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor v,
            FEVV::DataStructures::AIF::AIFMesh &sm)
{
  AIFHelpers::remove_edge(u, v, sm);
}

// BGL MutableEdgeGraph Concept [this function is needed when the data structure
// also follows BGL incidence graph concept]
inline void
remove_edge(typename boost::graph_traits<
                FEVV::DataStructures::AIF::AIFMesh >::out_edge_iterator
                iter, // mandatory to use out_edge_iterator (and not as we may
                      // think the usual edge_iterator)
            FEVV::DataStructures::AIF::AIFMesh &sm)
{
  sm.EraseIsolatedEdge(*iter);
}

// BGL MutableEdgeGraph Concept
template< typename UnaryPredicate >
inline void
remove_edge_if(UnaryPredicate p, FEVV::DataStructures::AIF::AIFMesh &sm)
{
  auto edges_range =
      (const_cast< FEVV::DataStructures::AIF::AIFMesh & >(sm)).GetEdges();
  std::vector< typename boost::graph_traits<
      FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor >
      edges_to_unlink;
  std::copy(
      edges_range.begin(),
      edges_range.end(),
      std::back_inserter(edges_to_unlink)); // the relations need to be copied
                                            // to keep iterator valid
  std::remove_if(edges_to_unlink.begin(), edges_to_unlink.end(), p);
}

// BGL MutableEdgeGraph Concept [this function is needed when the data structure
// also follows BGL incidence graph concept]
template< typename UnaryPredicate >
inline void
remove_out_edge_if(
    typename boost::graph_traits<
        FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor u,
    UnaryPredicate p,
    FEVV::DataStructures::AIF::AIFMesh &sm)
{
  auto edges_range_pair = out_edges(u, sm);
  std::vector< typename boost::graph_traits<
      FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor >
      edges_to_unlink;
  std::copy(
      edges_range_pair.first,
      edges_range_pair.second,
      std::back_inserter(edges_to_unlink)); // the relations need to be copied
                                            // to keep iterator valid
  std::remove_if(edges_to_unlink.begin(), edges_to_unlink.end(), p);
}

// BGL MutableEdgeGraph Concept [this function is needed when the data structure
// also follows BGL bidirectional graph concept]
template< typename UnaryPredicate >
inline void
remove_in_edge_if(typename boost::graph_traits<
                      FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor u,
                  UnaryPredicate p,
                  FEVV::DataStructures::AIF::AIFMesh &sm)
{
  auto edges_range_pair = in_edges(u, sm);
  std::vector< typename boost::graph_traits<
      FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor >
      edges_to_unlink;
  std::copy(
      edges_range_pair.first,
      edges_range_pair.second,
      std::back_inserter(edges_to_unlink)); // the relations need to be copied
                                            // to keep iterator valid
  std::remove_if(edges_to_unlink.begin(), edges_to_unlink.end(), p);
}

// BGL MutableEdgeGraph Concept
/// Remove all edges incident to v (to and from vertex v) from the mesh sm
/// Is usually called before remove_vertex(v, m) since this last method is not
/// supposed to do the unlinking with the incident edges
inline void
clear_vertex(typename boost::graph_traits<
                 FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor v,
             FEVV::DataStructures::AIF::AIFMesh &sm)
{
  AIFHelpers::clear_vertex(v);

  assert(AIFHelpers::is_isolated_vertex(v));
}

// CGAL MutableFaceGraph Concept
//!
//! \brief  Removes f from the mesh.
//!
inline void
remove_face(typename boost::graph_traits<
                FEVV::DataStructures::AIF::AIFMesh >::face_descriptor f,
            FEVV::DataStructures::AIF::AIFMesh &sm)
{
  AIFHelpers::remove_face(f, sm);
}

// CGAL MutableFaceGraph Concept
//!
//! \brief  Sets the corresponding face of h to f.
//!
inline void
set_face(typename boost::graph_traits< FEVV::DataStructures::AIF::AIFMesh >::
             halfedge_descriptor h, // cannot use a ref
         typename boost::graph_traits<
             FEVV::DataStructures::AIF::AIFMesh >::face_descriptor f,
         FEVV::DataStructures::AIF::AIFMesh &sm)
{
  AIFHelpers::set_face(h, f, sm);
}

// MUTS MutableFaceIncidentGraph concept refined by MutableCellIncidenceGraph
// Concept
inline void
add_in_edge(typename boost::graph_traits<
                FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor e,
            typename boost::graph_traits<
                FEVV::DataStructures::AIF::AIFMesh >::face_descriptor f)
{
  AIFHelpers::add_face(e, f);
}

// MUTS MutableFaceIncidentGraph concept refined by MutableCellIncidenceGraph
// Concept
inline void
remove_in_edge(typename boost::graph_traits<
                   FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor e,
               typename boost::graph_traits<
                   FEVV::DataStructures::AIF::AIFMesh >::face_descriptor f)
{
  AIFHelpers::remove_face(e, f);
}

// CGAL MutableFaceGraph Concept
//!
//! \brief  Sets the corresponding halfedge of f to h.
//!
inline void
set_halfedge(typename boost::graph_traits<
                 FEVV::DataStructures::AIF::AIFMesh >::face_descriptor f,
             typename boost::graph_traits<
                 FEVV::DataStructures::AIF::AIFMesh >::halfedge_descriptor h,
             FEVV::DataStructures::AIF::AIFMesh &sm)
{
  AIFHelpers::set_halfedge(f, h, sm);
}

// MUTS MutableFaceIncidentGraph concept refined by MutableCellIncidenceGraph
// Concept
inline void
add_out_edge(typename boost::graph_traits<
                 FEVV::DataStructures::AIF::AIFMesh >::face_descriptor f,
             typename boost::graph_traits<
                 FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor e)
{
  AIFHelpers::add_edge(f, e);
}

// MUTS MutableFaceIncidentGraph concept refined by MutableCellIncidenceGraph
// Concept
inline void
remove_out_edge(typename boost::graph_traits<
                    FEVV::DataStructures::AIF::AIFMesh >::face_descriptor f,
                typename boost::graph_traits<
                    FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor e)
{
  AIFHelpers::remove_edge(f, e);
}

// CGAL MutableFaceGraph Concept
//!
//! \brief  Adds a new face to the graph without initializing the connectivity.
//!
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::face_descriptor
add_face(FEVV::DataStructures::AIF::AIFMesh &sm)
{
  return AIFHelpers::add_face(sm);
}

} // namespace AIF
} // namespace DataStructures
} // namespace FEVV


// overload Euler::add_face(const VertexRange& vr, Graph& g) from
// Euler_operations.h ; this is not the most appropriate place to put this code
// because it is not related to any concept, but this is the most handy place to
// be sure that it will be found when Euler::add_face(vr, g) is used with AIF
namespace CGAL {
namespace Euler {
#ifndef USE_ADD_EDGE_AIF 
/**
* adds and returns the edge `e` connecting `s` and `t`
* halfedge(e, g) has s as source and t as target
*/
typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor
inline add_edge(typename boost::graph_traits<FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor s,
         typename boost::graph_traits<FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor t,
         FEVV::DataStructures::AIF::AIFMesh & g)
{
  typename boost::graph_traits<FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor e = AIFHelpers::add_edge(s, t, g).first;
  
  if (source(e, g) != s)
    AIFHelpers::swap_vertices(e);

  return e;
}
#endif
/**
* adds a new face defined by a range of vertices (identified by their descriptors,
* `boost::graph_traits<Graph>::%vertex_descriptor`).
* For each pair of consecutive vertices, the corresponding halfedge
* is added in `g` if new, and its connectivity is updated otherwise.
* The face can be added only at the boundary of `g`, or as a new connected component.
*
* @pre `vr` contains at least 3 vertices
* @returns the added face descriptor, or `boost::graph_traits<Graph>::%null_face()` if the face could not be added.
*/
template< typename VertexRange >
typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::face_descriptor
add_face(const VertexRange &vr, FEVV::DataStructures::AIF::AIFMesh &g)
{
  typedef FEVV::DataStructures::AIF::AIFMesh Graph;
  typedef typename boost::graph_traits< Graph >::vertex_descriptor
      vertex_descriptor;
  typedef
      typename boost::graph_traits< Graph >::face_descriptor face_descriptor;
  typedef
      typename boost::graph_traits< Graph >::edge_descriptor edge_descriptor;

  typedef FEVV::DataStructures::AIF::AIFTopologyHelpers AIFHelpers;

  std::vector< vertex_descriptor > vertices(vr.begin(),
                                            vr.end()); // quick and dirty copy
  unsigned int nv =
      (unsigned int)vertices.size(); // number of vertices in the face

  // build the face
  face_descriptor face = AIFHelpers::null_face();
#if 0
	//TODO-elo-choose  What to do if the face has 2 vertices only
	//                  -> solution 1: create one single edge and return a null face
	if( nv < 2 )
	{
		// don't allow degenerated faces
		throw std::runtime_error("add_face(vr, g) in Graph_traits_aif.h was called with less than 2 vertices");
	}
	else if( nv == 2 )
	{
		// Single edge

		// create the edge in the mesh
		edge_descriptor edge = AIFHelpers::add_edge(g);

		// link the edge and the vertices
		AIFHelpers::link_vertex_and_edge(vertices[0], edge, AIFHelpers::vertex_pos::FIRST);
		AIFHelpers::link_vertex_and_edge(vertices[1], edge, AIFHelpers::vertex_pos::SECOND);
	}
#else
  // TODO-elo-choose  What to do if the face has 2 vertices only
  //                  -> solution 2: do not create any edge, throw an exception
  //                  -> other solution (not implemented): do not create any
  //                  edge, return a null face
  if(nv < 3)
  {
    // don't allow degenerated faces
    throw std::runtime_error("add_face(vr, g) in Graph_traits_aif.h was called "
                             "with less than 3 vertices");
  }
#endif
  else if(nv > 2)
  {
    // create the face in the mesh
    face = AIFHelpers::add_face(g);

    // the source of the first edge is the last vertex in the list
    vertex_descriptor vsource = vertices.back();

    // loop over face's vertices to create the edges
    for(unsigned int i = 0; i < nv; i++)
    {
      vertex_descriptor vtarget = vertices[i];

      // create the edge if it doesn't already exist
      edge_descriptor edge = AIFHelpers::common_edge(vsource, vtarget);
      if(edge == AIFHelpers::null_edge())
      {
        // create the edge in the mesh
        edge = AIFHelpers::add_edge(g);

        // link the edge and the vertices
        AIFHelpers::link_vertex_and_edge(
            vsource, edge, AIFHelpers::vertex_pos::FIRST);
        AIFHelpers::link_vertex_and_edge(
            vtarget, edge, AIFHelpers::vertex_pos::SECOND);
      }

      // link the edge and the face (when needed -> else dangling edge)
	  if(! AIFHelpers::are_incident(edge, face) ||  ! AIFHelpers::are_incident(face, edge) )
        AIFHelpers::link_edge_and_face(edge, face);

      // the target of the current edge is the source of the next edge
      vsource = vtarget;
    }
  }

  return face;
}

} // namespace Euler
} // namespace CGAL


#if defined(BOOST_MSVC)
#pragma warning(pop)
#endif

