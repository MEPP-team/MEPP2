// Copyright (c) 2012-2022 University of Lyon and CNRS (France).
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
#include "Error_metric.h"

namespace FEVV {
namespace Filters {
template<
    typename HalfedgeGraph,
    typename PointMap>
class Volume_preserving : public Error_metric< HalfedgeGraph, PointMap >
{
public:
  using edge_iterator =
            typename boost::graph_traits< HalfedgeGraph >::edge_iterator;
  using edge_descriptor =
            typename boost::graph_traits< HalfedgeGraph >::edge_descriptor;
  using Vector = typename FEVV::Geometry_traits< HalfedgeGraph >::Vector;
  using Point = typename FEVV::Geometry_traits< HalfedgeGraph >::Point;
  using Geometry = typename FEVV::Geometry_traits< HalfedgeGraph>;
  
  typedef Error_metric< HalfedgeGraph, PointMap > Super_class;
  
  Volume_preserving(HalfedgeGraph &g,
                   PointMap &pm,
                   Kept_position< HalfedgeGraph, PointMap > *vkept,
      FEVV::Filters::Uniform_dequantization< HalfedgeGraph, PointMap > &dequantiz)
      : Super_class(g, pm, vkept, dequantiz){}
	  
  ~Volume_preserving(){}
  
  void compute_error() override
  {
    if(!Super_class::_queue.empty())
    {
      typename Super_class::priority_queue_edges empty2;
      std::swap(Super_class::_queue, empty2);
    }
    Super_class::_edges_cost.clear();
    typename Super_class::edge2cost_map empty;
    std::swap(Super_class::_edges_cost, empty);
    if(Super_class::_edges_cost.empty())
    {
      //std::cout << "Volume Metric" << std::endl;
      //Super_class::_pm = get(boost::vertex_point, Super_class::_g);
      auto edge_iterator_pair = edges(Super_class::_g);
      auto edge_ite = edge_iterator_pair.first;
      int count = 0;
      Super_class::_threshold = 0;
      for( ; edge_ite != edge_iterator_pair.second;
          ++edge_ite) 
      {
        Point collapsePos = Super_class::_vkept->compute_position(*edge_ite);

        double weight = compute_cost_edge(*edge_ite, collapsePos);
        Super_class::_threshold += weight;
        ++count;
		
        Super_class::_queue.push(
                  std::make_tuple(*edge_ite, weight, collapsePos));
        Super_class::_edges_cost.emplace(*edge_ite,
                             std::make_pair(weight, collapsePos));
      }
      Super_class::_threshold /= count;
    }
  };

  double compute_cost_edge(edge_descriptor e) 
  {
    Point collapsePos = Super_class::_vkept->compute_position(e);
	  
    double sum_first_vertex=0., sum_second_vertex=0.;
    compute_volumes(e, collapsePos, sum_first_vertex, sum_second_vertex);
	
    return (sum_first_vertex + sum_second_vertex);
  }

  double compute_cost_edge(edge_descriptor e, const Point &collapsePos) override
  {
    double sum_first_vertex=0., sum_second_vertex=0.;
    compute_volumes(e, collapsePos, sum_first_vertex, sum_second_vertex);
    
    return (sum_first_vertex + sum_second_vertex);
  }
  std::string get_as_string() const override
  {
	  return "volumepreserving";
  }
  private:
  double tetrahedron_volume(const Point &a,
                           const Point &b,
                           const Point &c,
                           const Point &d) const 
  {
    return fabs(Super_class::_gt.dot_product(
               Super_class::_gt.cross_product(Super_class::_gt.sub_p(b, a), Super_class::_gt.sub_p(c, a)),
               Super_class::_gt.sub_p(d, a))) /
           6.0; // this calculus has been validated onto several
    // tetrahedra using the above commented source code
  };
  void compute_volumes(edge_descriptor e, const Point &collapsePos, double& sum_first_vertex, double& sum_second_vertex)
  {
    auto current_halfedge = halfedge(e, Super_class::_g);
    boost::iterator_range<
        CGAL::Halfedge_around_target_iterator< HalfedgeGraph > >
        iterator_range = CGAL::halfedges_around_target(target(current_halfedge, Super_class::_g), Super_class::_g);
    sum_first_vertex = 0.;
    for(auto h : iterator_range)
    {	  
      if((h == current_halfedge) ||
         (prev(opposite(h, Super_class::_g), Super_class::_g) == current_halfedge) ||
         CGAL::is_border(opposite(h, Super_class::_g),
                         Super_class::_g)) // when volume calculus is not possible
                                          // we ignore it
      {
        continue;
      }
      const Point& a = get(Super_class::_pm, source(h, Super_class::_g));
      Point b = get(Super_class::_pm, target(next(opposite(h, Super_class::_g), Super_class::_g), Super_class::_g));
      const Point& c = get(Super_class::_pm, target(current_halfedge, Super_class::_g));

      sum_first_vertex +=
          tetrahedron_volume(a, b, c, collapsePos); // volume polyhedron
    }
    iterator_range =
        CGAL::halfedges_around_target(source(current_halfedge, Super_class::_g), Super_class::_g);
    sum_second_vertex = 0.;
    for(auto h : iterator_range)
    {
      if((h == opposite(current_halfedge, Super_class::_g)) ||
         (prev(opposite(h, Super_class::_g), Super_class::_g) == opposite(current_halfedge, Super_class::_g)) ||
         (CGAL::is_border(opposite(h, Super_class::_g), Super_class::_g)))
      {
        continue;
      }
      const Point& a = get(Super_class::_pm, source(h, Super_class::_g));
      Point b = get(Super_class::_pm, target(next(opposite(h, Super_class::_g), Super_class::_g), Super_class::_g));
      const Point& c = get(Super_class::_pm, source(current_halfedge, Super_class::_g));

      sum_second_vertex +=
          tetrahedron_volume(a,
                            b,
                            c,
                            collapsePos);
    }
  }
};
} // namespace Filters
} // namespace FEVV