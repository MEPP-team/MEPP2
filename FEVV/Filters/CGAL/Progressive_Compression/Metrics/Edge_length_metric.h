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
#include <iostream>

namespace FEVV {
namespace Filters {
template<
	typename HalfedgeGraph,
	typename PointMap >
class Edge_length_metric : public Error_metric< HalfedgeGraph, PointMap >
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
  Edge_length_metric(HalfedgeGraph &g,
                   PointMap &pm,
                   Kept_position< HalfedgeGraph, PointMap > *vkept,
      FEVV::Filters::Uniform_dequantization< HalfedgeGraph, PointMap > &dequantiz)
      : Super_class(g, pm, vkept, dequantiz)
  {
    Super_class::_operator = FEVV::Filters::VKEPT_POSITION::MIDPOINT;
  }
  
  ~Edge_length_metric(){}
  
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
      //std::cout << "Edge length Metric" << std::endl;
      auto edge_iterator_pair = edges(Super_class::_g);
      auto edge_ite = edge_iterator_pair.first;
      int count = 0;
      Super_class::_threshold = 0;
      for( ; edge_ite != edge_iterator_pair.second; ++edge_ite)
      {
        Point collapsePos = Super_class::_vkept->compute_position(*edge_ite);

        double weight = compute_cost_edge(*edge_ite,collapsePos);

        Super_class::_threshold += weight;
        count++;

        Super_class::_queue.push(
                  std::make_tuple(*edge_ite, weight, collapsePos));
        Super_class::_edges_cost.emplace(*edge_ite,
                               std::make_pair(weight,collapsePos));
      }

     Super_class::_threshold /= count;
    }
  }

  std::string get_as_string() const override { return "edgelength";}

  double compute_cost_edge(edge_descriptor e)
  {
    const Point& P0 = get(Super_class::_pm, source(e, Super_class::_g));
    const Point& P1 = get(Super_class::_pm, target(e, Super_class::_g));
    return Super_class::_gt.length(Super_class::_gt.sub_p(P1, P0));
  }

  double compute_cost_edge(edge_descriptor e, const Point &/*collapsePos*/) override
  {
    return compute_cost_edge(e);
  }
};


} // namespace Filters
} // namespace FEVV
