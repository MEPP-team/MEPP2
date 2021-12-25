#pragma once

#include "ErrorMetric.h"
#include <iostream>

namespace FEVV {
namespace Filters {
template<
	typename HalfedgeGraph,
	typename PointMap >
class EdgeLengthMetric : public ErrorMetric< HalfedgeGraph, PointMap >
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
  EdgeLengthMetric(HalfedgeGraph &g,
                   PointMap &pm,
                   KeptPosition< HalfedgeGraph, PointMap > *vkept,
      FEVV::Filters::UniformDequantization< HalfedgeGraph, PointMap > &dequantiz)
      : SuperClass(g, pm, vkept, dequantiz)
  {
    SuperClass::_operator = FEVV::Filters::VKEPT_POSITION::MIDPOINT;
  }
  
  ~EdgeLengthMetric(){}
  
  void ComputeError() override
  {
    if(!SuperClass::_queue.empty())
    {
      typename SuperClass::priority_queue_edges empty2;
      std::swap(SuperClass::_queue, empty2);
    }
    SuperClass::_edges_cost.clear();
    if(SuperClass::_edges_cost.empty())
    {
      //std::cout << "Edge length Metric" << std::endl;
      auto edge_iterator_pair = edges(SuperClass::_g);
      auto edge_ite = edge_iterator_pair.first;
      int count = 0;
      SuperClass::_threshold = 0;
      for(; edge_ite != edge_iterator_pair.second; ++edge_ite)
      {
        Point collapsePos = SuperClass::_vkept->ComputePosition(*edge_ite);

        double weight = ComputeCostEdge(*edge_ite,collapsePos);

        SuperClass::_threshold += weight;
        count++;

        SuperClass::_queue.push(
                  std::make_tuple(*edge_ite, weight, collapsePos));
        SuperClass::_edges_cost.emplace(*edge_ite,
                               std::make_pair(weight,collapsePos));
      }

     SuperClass::_threshold /= count;
    }
  }

  std::string getMethodasString() const override { return "edgelength";}

  double ComputeCostEdge(edge_descriptor e)
  {
    const Point& P0 = get(SuperClass::_pm, source(e, SuperClass::_g));
    const Point& P1 = get(SuperClass::_pm, target(e, SuperClass::_g));
    return SuperClass::_gt.length(SuperClass::_gt.sub_p(P1, P0));
  }

  double ComputeCostEdge(edge_descriptor e, const Point &/*collapsePos*/) override
  {
    return ComputeCostEdge(e);
  }
};


} // namespace Filters
} // namespace FEVV
