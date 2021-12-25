#pragma once

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <CGAL/boost/graph/iterator.h>

namespace FEVV {
namespace Filters {
template< typename HalfedgeGraph,
          typename vertex_descriptor =
              typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor,
          typename halfedge_descriptor = typename boost::graph_traits<
              HalfedgeGraph >::halfedge_descriptor >

void
ForbidVertex(HalfedgeGraph &g,
             vertex_descriptor v,
             std::set< halfedge_descriptor > &forbidden_edges,
             std::vector< halfedge_descriptor > &/*edges_to_color*/)
{
  boost::iterator_range<
      CGAL::Halfedge_around_target_iterator< HalfedgeGraph > >
      iterator_range = CGAL::halfedges_around_target(v, g);
  for(auto h_v : iterator_range)
  {
    forbidden_edges.insert(h_v);
    forbidden_edges.insert(opposite(h_v, g));
  }
}

template< typename HalfedgeGraph,
          typename vertex_descriptor =
              typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor,
          typename halfedge_descriptor = typename boost::graph_traits<
              HalfedgeGraph >::halfedge_descriptor >
void
FindVerticesToForbid(HalfedgeGraph &g,
                     vertex_descriptor v,
                     std::set< halfedge_descriptor > &/*forbidden_edges*/,
                     std::set< vertex_descriptor > &forbidden_vertices,
                     std::vector< halfedge_descriptor > &/*edges_to_color*/)
{
  boost::iterator_range<
      CGAL::Halfedge_around_target_iterator< HalfedgeGraph > >
      iterator_range = CGAL::halfedges_around_target(v, g);
  forbidden_vertices.insert(v);
  for(auto h_v : iterator_range)
  {
    vertex_descriptor sourc_h = source(h_v, g);
    forbidden_vertices.insert(sourc_h);
    if(!CGAL::is_border(h_v, g))
    {

      halfedge_descriptor next_h = next(opposite(h_v, g), g);
      halfedge_descriptor opp_next_h = opposite(next_h, g);
      if(!CGAL::is_border(opp_next_h, g))
      {

        //halfedge_descriptor wing_h = next(opp_next_h, g);
        //forbidden_vertices.insert(target(wing_h, g));
      }
    }
  }
}


template< typename HalfedgeGraph,
          typename vertex_descriptor =
              typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor,
          typename halfedge_descriptor = typename boost::graph_traits<
              HalfedgeGraph >::halfedge_descriptor >
void
ForbidEdges(HalfedgeGraph &g,
            vertex_descriptor v,
            std::set< halfedge_descriptor > &forbidden_edges,
            std::vector< halfedge_descriptor > &edges_to_color)
{
  std::set< vertex_descriptor > forbidden_vertices;
  FindVerticesToForbid(
      g, v, forbidden_edges, forbidden_vertices, edges_to_color);
  for(const vertex_descriptor &ver : forbidden_vertices)
  {
    ForbidVertex(g, ver, forbidden_edges, edges_to_color);
  }
}

template< typename HalfedgeGraph,
          typename vertex_descriptor =
              typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor,
          typename halfedge_descriptor = typename boost::graph_traits<
              HalfedgeGraph >::halfedge_descriptor >
void
ComputeSimpleStencil(HalfedgeGraph &g,
                     vertex_descriptor v,
                     std::set< halfedge_descriptor > &forbidden_edges,
                     std::vector< halfedge_descriptor > &edges_to_color)
{
  boost::iterator_range<
      CGAL::Halfedge_around_target_iterator< HalfedgeGraph > >
      iterator_range = CGAL::halfedges_around_target(v, g);
  for(auto h_v : iterator_range)
  {
    forbidden_edges.insert(h_v);
    forbidden_edges.insert(opposite(h_v, g));

    edges_to_color.push_back(h_v);
    edges_to_color.push_back(opposite(h_v, g));

    halfedge_descriptor next_h = next(opposite(h_v, g), g);

    forbidden_edges.insert(next_h);
    halfedge_descriptor opp_next_h = opposite(next_h, g);
    forbidden_edges.insert(opp_next_h);
    edges_to_color.push_back(next_h);
    edges_to_color.push_back(opp_next_h);


    if(!CGAL::is_border(opp_next_h, g) &&
       !CGAL::is_border(next_h, g)) // second verification might be overkill
    {


      forbidden_edges.insert(next(opp_next_h, g));
      vertex_descriptor vv = target(opp_next_h, g);
      boost::iterator_range<
          CGAL::Halfedge_around_target_iterator< HalfedgeGraph > >
          iterator_range2 = CGAL::halfedges_around_target(vv, g);
      for(auto h_v2 : iterator_range2)
      {
        forbidden_edges.insert(h_v2);
        forbidden_edges.insert(opposite(h_v2, g));
      }
      forbidden_edges.insert(opposite(next(opp_next_h, g), g));

      edges_to_color.push_back(next(opp_next_h, g));
      edges_to_color.push_back(opposite(next(opp_next_h, g), g));

      forbidden_edges.insert(next(next(opp_next_h, g), g));
      forbidden_edges.insert(opposite(next(next(opp_next_h, g), g), g));

      edges_to_color.push_back(next(next(opp_next_h, g), g));
      edges_to_color.push_back(opposite(next(next(opp_next_h, g), g), g));
    }
  }
}

} // namespace Filters
} // namespace FEVV
