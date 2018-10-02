#pragma once

#if defined(SplitEdgeNonManifoldFilter_RECURSES)
#error Recursive header files inclusion detected in SplitEdgeNonManifoldFilter.h
#else // defined(SplitEdgeNonManifoldFilter_RECURSES)
/** Prevents recursive inclusion of headers. */
#define SplitEdgeNonManifoldFilter_RECURSES

#if !defined SplitEdgeNonManifoldFilter_h
/** Prevents repeated inclusion of headers. */
#define SplitEdgeNonManifoldFilter_h
//#include <CGAL/boost/graph/internal/helpers.h> // set_border
#include <boost/graph/graph_traits.hpp>
#include <CGAL/boost/graph/iterator.h>

#include <CGAL/boost/graph/Euler_operations.h> // for add_face(vr,g)
//#include <CGAL/boost/graph/helpers.h> // is_border + is_triangle_mesh

#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Operators/AIF/topology_predicates.hpp"

#include <cassert>

namespace FEVV {
namespace Filters {

/**
 * \brief Split an edge of the graph.
 *
 * \tparam  MutableFaceIncidentGraph a Mesh type that provides a Model of the
 *          MutableFaceIncidentGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \tparam  PointMap A modifiable point map to manage vertex positions.
 * \param g The MutableFaceIncidentGraph instance from which the e edge will be
 *          split.
 * \param pm The mesh point map which associates vertex to positions.
 * \param e The edge to be split.
 */
template< typename MutableFaceIncidentGraph, // similar concept to
                                             // MutableCellIncidentGraph, but
                                             // limited to 0, 1 and 2
                                             // dimensional cells
          typename PointMap >
void
split_edge(
    MutableFaceIncidentGraph &g,
    PointMap pm,
    typename boost::graph_traits< MutableFaceIncidentGraph >::edge_descriptor
        &e)
{
  typedef FEVV::Geometry_traits< MutableFaceIncidentGraph > GeometryTraits;
  typedef typename GeometryTraits::Point Point;

  typedef boost::graph_traits< MutableFaceIncidentGraph > GraphTraits;
  typedef typename GraphTraits::edge_descriptor edge_descriptor;
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
  typedef typename GraphTraits::face_descriptor face_descriptor;
  GeometryTraits gt(g);

  if(e == boost::graph_traits< MutableFaceIncidentGraph >::null_edge())
    return;

  // ensure the edge is incident to triangles (only)
  if(!FEVV::Operators::has_only_incident_triangular_faces(e, g))
  {
    return;
  }

  // deal with vertex geometry
  vertex_descriptor vs = source(e, g);
  vertex_descriptor vt = target(e, g);

  vertex_descriptor midpoint_vertex;
  midpoint_vertex = add_vertex(g);
  put(pm,
      midpoint_vertex,
      Point((gt.get_x(get(pm, vs)) + gt.get_x(get(pm, vt))) * 0.5f,
            (gt.get_y(get(pm, vs)) + gt.get_y(get(pm, vt))) * 0.5f,
            (gt.get_z(get(pm, vs)) + gt.get_z(get(pm, vt))) * 0.5f));

  // when a topological modication is done over an edge, better to remove that
  // edge than to reuse it because we have to make sure the final user will
  // update edge properties accordingly.
  if(degree(e, g) == 0)
  { // one dangling edge is tranfromed into 2 successive dangling edges
    add_edge(vs, midpoint_vertex, g);
    add_edge(midpoint_vertex, vt, g);

    remove_edge(e, g); // remove at the end, to make sure no vertex is removed
  }
  else
  {
    std::vector< vertex_descriptor > face1_vertices, face2_vertices;
    auto face_range_pair = in_edges(e, g); // get incident faces
    std::vector< face_descriptor > copy_f(
        face_range_pair.first,
        face_range_pair
            .second); // copy to not invalidate iterator while iterating
    typename std::vector< face_descriptor >::iterator iter_f(copy_f.begin()),
        iter_fe(copy_f.end());
    for(; iter_f != iter_fe; ++iter_f)
    {
      auto edge_range_pair = out_edges(*iter_f, g);
      auto iter_e = edge_range_pair.first;
      vertex_descriptor third_v;
      /////////////////////////////////////////////////////////////////////////
      while(*iter_e != e)
        ++iter_e;
      /////////////////////////////////////////////////////////////////////////
      edge_descriptor prev_e = *iter_e;
      ++iter_e;
      if(iter_e == edge_range_pair.second)
        iter_e = edge_range_pair.first;
      /////////////////////////////////////////////////////////////////////////
      if(source(prev_e, g) == source(*iter_e, g))
      {
        face1_vertices.push_back(target(*iter_e, g));
        face1_vertices.push_back(target(prev_e, g));
        face1_vertices.push_back(midpoint_vertex);

        face2_vertices.push_back(source(prev_e, g));
        face2_vertices.push_back(target(*iter_e, g));
        face2_vertices.push_back(midpoint_vertex);
      }
      else if(target(prev_e, g) == source(*iter_e, g))
      {
        face1_vertices.push_back(target(*iter_e, g));
        face1_vertices.push_back(source(prev_e, g));
        face1_vertices.push_back(midpoint_vertex);

        face2_vertices.push_back(target(prev_e, g));
        face2_vertices.push_back(target(*iter_e, g));
        face2_vertices.push_back(midpoint_vertex);
      }
      else if(source(prev_e, g) == target(*iter_e, g))
      {
        face1_vertices.push_back(source(*iter_e, g));
        face1_vertices.push_back(target(prev_e, g));
        face1_vertices.push_back(midpoint_vertex);

        face2_vertices.push_back(source(prev_e, g));
        face2_vertices.push_back(source(*iter_e, g));
        face2_vertices.push_back(midpoint_vertex);
      }
      else if(target(prev_e, g) == target(*iter_e, g))
      {
        face1_vertices.push_back(source(*iter_e, g));
        face1_vertices.push_back(source(prev_e, g));
        face1_vertices.push_back(midpoint_vertex);

        face2_vertices.push_back(target(prev_e, g));
        face2_vertices.push_back(source(*iter_e, g));
        face2_vertices.push_back(midpoint_vertex);
      }

      CGAL::Euler::add_face(face1_vertices, g);
      CGAL::Euler::add_face(face2_vertices, g);
      face1_vertices.clear();
      face2_vertices.clear();

      remove_face(*iter_f,
                  g); // remove at the end, to make sure no vertex or edge is
                      // removed while still needing by new geometries
    }
  }
}

} // namespace Filters
} // namespace FEVV

#endif // !defined SplitEdgeNonManifoldFilter_h

#undef SplitEdgeNonManifoldFilter_RECURSES
#endif // else defined(SplitEdgeNonManifoldFilter_RECURSES)
