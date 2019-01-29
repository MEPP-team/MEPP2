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

#include <map>
#include <vector>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <CGAL/boost/graph/helpers.h> // is_border
#include "FEVV/Wrappings/Geometry_traits.h"

namespace FEVV {
namespace Operators {

 /**
 * \brief Extract vertices which are at most k (inclusive)
 *        far from vertex v in the graph of edges.
 * \tparam  FaceGraph a Mesh type that provides a Model of the
 *          FaceGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \tparam  GeometryTraits The geometric kernel when available. This is
 *          defaulted to FEVV::Geometry_traits<FaceGraph>.
 * \param[in]  v The vertex whose k-rings are computed. 
 * \param[in]  g The FaceGraph instance.
 * \param[in]  k The maximum distance to include a vertex v' neighbor. 
 * \param[out] qv The std::vector of vertices at most k (inclusive)
 *             far from vertex v in the graph of edges
 */ 
template< typename FaceGraph,
          typename GeometryTraits = FEVV::Geometry_traits< FaceGraph > >
void
extract_k_ring(
    typename boost::graph_traits< FaceGraph >::vertex_descriptor v,
    const FaceGraph &g,
    int k,
    std::vector< typename boost::graph_traits< FaceGraph >::vertex_descriptor >
        &qv)
{
  typedef boost::graph_traits< FaceGraph > GraphTraits;
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
  typedef typename GraphTraits::face_descriptor face_descriptor;

  std::map< vertex_descriptor, int > d;
  qv.push_back(v);
  d[v] = 0;
  std::size_t current_index = 0;
  int dist_v;
  while(current_index < qv.size() && (dist_v = d[qv[current_index]]) < k)
  {
    v = qv[current_index++];
    halfedge_descriptor h = halfedge(v, g), hend = h;
    // CGAL::Vertex_around_target_circulator<HalfedgeGraph> e(halfedge(v, g),
    // g), e_end(e); // for surface mesh
    // HalfedgeGraph::Halfedge_around_vertex_circulator e(v->vertex_begin()),
    // e_end(e); // ok for polyhedron
    do
    {
      vertex_descriptor new_v = source(h, g), tardeb = target(h, g);
      assert(v == tardeb);
      if(d.insert(std::make_pair(new_v, dist_v + 1)).second)
      {
        qv.push_back(new_v);
      }
      // There is a pb when the next edge is a border edge: we must take into
      // account the vertices while updating correctly the next halfedge
      if(face(opposite(next(h, g), g), g) ==
         boost::graph_traits< FaceGraph >::null_face())
      {
        if(d.insert(std::make_pair(target(next(h, g), g), dist_v + 1)).second)
        {
          qv.push_back(target(next(h, g), g));
        }
      }
      h = opposite(next(h, g), g);
    } while(!(h == hend));
  }
  // debug
  // std::cout << "qv.size() = " << qv.size() << std::endl;
}


 /**
 * \brief Extract halfedges which have v as target vertex.
 *        All extracted halfedges are associated to a face, 
 *        meaning that for a border halfedge, it is its 
 *        opposite halfedge that is extracted.
 * \tparam  FaceGraph a Mesh type that provides a Model of the
 *          FaceGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \param[in]  v The vertex whose incident halfedges are extracted. 
 * \param[in]  g The FaceGraph instance.
 * \param[out] qh The std::vector of halfedges incident to v.
 */  
template< typename FaceGraph >
void
extract_vertex_star(
    typename boost::graph_traits< FaceGraph >::vertex_descriptor v,
    const FaceGraph &g,
    std::vector<
        typename boost::graph_traits< FaceGraph >::halfedge_descriptor > &qh)
{
  typedef boost::graph_traits< FaceGraph > GraphTraits;
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
  typedef typename GraphTraits::face_descriptor face_descriptor;

  std::size_t current_index = 0;

  halfedge_descriptor h = halfedge(v, g), hend = h;
  do
  {
    // vertex_descriptor new_v = source(h, g), tardeb = target(h, g);
    // assert(v == tardeb);
    qh.push_back(h);
    vertex_descriptor vs = source(h, g), vt = target(h, g); // debug
    // There is a pb when the next edge is a border edge: we must take into
    // account the vertices while updating correctly the next halfedge
    if(face(opposite(next(h, g), g), g) ==
       boost::graph_traits< FaceGraph >::null_face())
    {
      if(opposite(next(h, g), g) == qh[0]) // Surface Mesh
      {
        qh.erase(qh.begin());
        qh.push_back(opposite(halfedge(v, g), g));
        break;
      }
      qh.push_back(next(h, g)); // border edges are put in the other direction
                                // (because the opposite halfedge was null for
                                // AIFMesh first release [no more the case])
      // vs = source(next(h, g), g); vt = target(next(h, g), g); // debug
      // halfedge_descriptor htmp = opposite(next(h, g), g);
      if(!CGAL::is_border(
             opposite(h, g),
             g) // "not a border edge" for polyhedron, OpenMesh and AIF

         //&& !(opposite(h, g) ==
         //boost::graph_traits<FaceGraph>::null_halfedge()) // no condition
         //works properly with Surface Mesh
         //! CGAL::is_border(target(opposite(h, g), g),g)
         ) // manage the case where there is only one valid halfedge
      {
        halfedge_descriptor last = h;
        // we try to reach the other border (simple, because we assume manifold
        // configuration)
        do
        {
          if(prev(opposite(h, g), g) == last) // for surface mesh...
            break;
          h = prev(opposite(h, g), g);
          // vs = source(h, g); vt = target(h, g); // debug
        } while(!CGAL::is_border(
            opposite(h, g),
            g) // works with polyhedron, OpenMesh and AIF
               // !(opposite(h, g) ==
               // boost::graph_traits<FaceGraph>::null_halfedge()) // no
               // condition works properly with Surface Mesh
        );
      }
    }
    else
      h = opposite(next(h, g), g);
  } while(!(h == hend));
  // debug
  // std::cout << "qh.size() = " << qh.size() << std::endl;
}

 /**
 * \brief Extract vertices which are at exactly 1
 *        far from vertex v in the graph of edges.
 * \tparam  FaceGraph a Mesh type that provides a Model of the
 *          FaceGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \tparam  GeometryTraits The geometric kernel when available. This is
 *          defaulted to FEVV::Geometry_traits<FaceGraph>.
 * \param[in]  v The vertex whose 1-ring is computed. 
 * \param[in]  g The FaceGraph instance.
 * \param[out] qv The std::vector of vertices at exactly 1
 *             far from vertex v in the graph of edges
 */ 
template< typename FaceGraph,
          typename GeometryTraits = FEVV::Geometry_traits< FaceGraph > >
void
extract_1_ring_not_including_v(
    typename boost::graph_traits< FaceGraph >::vertex_descriptor v,
    const FaceGraph &g,
    std::vector< typename boost::graph_traits< FaceGraph >::vertex_descriptor >
        &qv)
{
  Operators::extract_k_ring< FaceGraph, GeometryTraits >(v, g, 1, qv);
  qv.erase(qv.begin()); // the first element is v itself ant it must be removed
}

} // namespace Operators
} // namespace FEVV
