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
#include "ErrorMetric.h"

namespace FEVV {
namespace Filters {
template<
    typename HalfedgeGraph,
    typename PointMap>
class VolumePreserving : public ErrorMetric< HalfedgeGraph, PointMap >
{
public:
  using edge_iterator =
            typename boost::graph_traits< HalfedgeGraph >::edge_iterator;
  using edge_descriptor =
            typename boost::graph_traits< HalfedgeGraph >::edge_descriptor;
  using Vector = typename FEVV::Geometry_traits< HalfedgeGraph >::Vector;
  using Point = typename FEVV::Geometry_traits< HalfedgeGraph >::Point;
  using Geometry = typename FEVV::Geometry_traits< HalfedgeGraph>;
  
  typedef ErrorMetric< HalfedgeGraph, PointMap > SuperClass;
  
  VolumePreserving(HalfedgeGraph &g,
                   PointMap &pm,
                   KeptPosition< HalfedgeGraph, PointMap > *vkept,
      FEVV::Filters::UniformDequantization< HalfedgeGraph, PointMap > &dequantiz)
      : SuperClass(g, pm, vkept, dequantiz){}
	  
  ~VolumePreserving(){}
  
  void ComputeError() override
  {
    if (!SuperClass::_queue.empty())
    {
      typename SuperClass::priority_queue_edges empty2;
      std::swap(SuperClass::_queue, empty2);
    }
    SuperClass::_edges_cost.clear();
    typename SuperClass::edge2cost_map empty;
    std::swap(SuperClass::_edges_cost, empty);
    if(SuperClass::_edges_cost.empty())
    {
      //std::cout << "Volume Metric" << std::endl;
      //SuperClass::_pm = get(boost::vertex_point, SuperClass::_g);
      auto edge_iterator_pair = edges(SuperClass::_g);
      auto edge_ite = edge_iterator_pair.first;
      int count = 0;
      SuperClass::_threshold = 0;
      for(; edge_ite != edge_iterator_pair.second;
          ++edge_ite) 
      {
        Point collapsePos = SuperClass::_vkept->ComputePosition(*edge_ite);

        double weight = ComputeCostEdge(*edge_ite, collapsePos);
        SuperClass::_threshold += weight;
        ++count;
		
        SuperClass::_queue.push(
                  std::make_tuple(*edge_ite, weight, collapsePos));
        SuperClass::_edges_cost.emplace(*edge_ite,
                             std::make_pair(weight, collapsePos));
      }
      SuperClass::_threshold /= count;
    }
  };

  double ComputeCostEdge(edge_descriptor e) 
  {
    Point collapsePos = SuperClass::_vkept->ComputePosition(e);
	  
    double sum_first_vertex=0., sum_second_vertex=0.;
    compute_volumes(e, collapsePos, sum_first_vertex, sum_second_vertex);
	
    return (sum_first_vertex + sum_second_vertex);
  }

  double ComputeCostEdge(edge_descriptor e, const Point &collapsePos) override
  {
    double sum_first_vertex=0., sum_second_vertex=0.;
    compute_volumes(e, collapsePos, sum_first_vertex, sum_second_vertex);
    
    return (sum_first_vertex + sum_second_vertex);
  }
  std::string getMethodasString() const override
  {
	  return "volumepreserving";
  }
  private:
  double tetrahedronVolume(const Point &a,
                           const Point &b,
                           const Point &c,
                           const Point &d) const 
  {
    return fabs(SuperClass::_gt.dot_product(
               SuperClass::_gt.cross_product(SuperClass::_gt.sub_p(b, a), SuperClass::_gt.sub_p(c, a)),
               SuperClass::_gt.sub_p(d, a))) /
           6.0; // this calculus has been validated onto several
    // tetrahedra using the above commented source code
  };
  void compute_volumes(edge_descriptor e, const Point &collapsePos, double& sum_first_vertex, double& sum_second_vertex)
  {
    auto current_halfedge = halfedge(e, SuperClass::_g);
    boost::iterator_range<
        CGAL::Halfedge_around_target_iterator< HalfedgeGraph > >
        iterator_range = CGAL::halfedges_around_target(target(current_halfedge, SuperClass::_g), SuperClass::_g);
    sum_first_vertex = 0.;
    for(auto h : iterator_range)
    {	  
      if((h == current_halfedge) ||
         (prev(opposite(h, SuperClass::_g), SuperClass::_g) == current_halfedge) ||
         CGAL::is_border(opposite(h, SuperClass::_g),
                         SuperClass::_g)) // when volume calculus is not possible
                                          // we ignore it
      {
        continue;
      }
      const Point& a = get(SuperClass::_pm, source(h, SuperClass::_g));
      Point b = get(SuperClass::_pm, target(next(opposite(h, SuperClass::_g), SuperClass::_g), SuperClass::_g));
      const Point& c = get(SuperClass::_pm, target(current_halfedge, SuperClass::_g));

      sum_first_vertex +=
          tetrahedronVolume(a, b, c, collapsePos); // volume polyhedron
    }
    iterator_range =
        CGAL::halfedges_around_target(source(current_halfedge, SuperClass::_g), SuperClass::_g);
    sum_second_vertex = 0.;
    for(auto h : iterator_range)
    {
      if((h == opposite(current_halfedge, SuperClass::_g)) ||
         (prev(opposite(h, SuperClass::_g), SuperClass::_g) == opposite(current_halfedge, SuperClass::_g)) ||
         (CGAL::is_border(opposite(h, SuperClass::_g), SuperClass::_g)))
      {
        continue;
      }
      const Point& a = get(SuperClass::_pm, source(h, SuperClass::_g));
      Point b = get(SuperClass::_pm, target(next(opposite(h, SuperClass::_g), SuperClass::_g), SuperClass::_g));
      const Point& c = get(SuperClass::_pm, source(current_halfedge, SuperClass::_g));

      sum_second_vertex +=
          tetrahedronVolume(a,
                            b,
                            c,
                            collapsePos);
    }
  }
};
} // namespace Filters
} // namespace FEVV