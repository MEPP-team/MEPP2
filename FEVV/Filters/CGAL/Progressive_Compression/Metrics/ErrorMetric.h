#pragma once

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <CGAL/boost/graph/iterator.h>
#include "FEVV/Wrappings/Geometry_traits.h"

#include "FEVV/Filters/CGAL/Progressive_Compression/Operators/Midpoint.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Operators/Halfedge.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/dequantization.h"

#include <iostream>
#include <vector>
#include <map>
#include <queue>
#include <utility>

namespace FEVV {
namespace Filters {
template<
	typename HalfedgeGraph,
	typename edge_descriptor =
	typename boost::graph_traits< HalfedgeGraph >::edge_descriptor,
	typename Point = typename FEVV::Geometry_traits< HalfedgeGraph >::Point,
	typename Vector = typename FEVV::Geometry_traits< HalfedgeGraph >::Vector >
class Compare_weights2
{
public:
  bool operator()(const std::tuple< edge_descriptor, double, Point > &e1,
                  const std::tuple< edge_descriptor, double, Point > &e2) const
  {
    double w1 = std::get< 1 >(e1);
    double w2 = std::get< 1 >(e2);

    if(w1 > w2)
      return true;
    else
      return false;
  }
};

/// \brief Abstract class to compute the collapse cost of each edge in a mesh.
template<
    typename HalfedgeGraph,
    typename PointMap>
class ErrorMetric
{
public:
  using vertex_descriptor =
      typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor;
  using halfedge_descriptor =
      typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor;
        using edge_iterator =
            typename boost::graph_traits< HalfedgeGraph >::edge_iterator;
  using edge_descriptor =
            typename boost::graph_traits< HalfedgeGraph >::edge_descriptor;
  using Vector = typename FEVV::Geometry_traits< HalfedgeGraph >::Vector;
  using Point = typename FEVV::Geometry_traits< HalfedgeGraph >::Point;
  using Geometry = typename FEVV::Geometry_traits< HalfedgeGraph>;


  typedef std::priority_queue<
      std::tuple< edge_descriptor, double, Point >,
      std::vector< std::tuple< edge_descriptor, double, Point > >, 
      Compare_weights2< HalfedgeGraph > >
    priority_queue_edges;

  typedef std::map< edge_descriptor, std::pair< double, Point > > edge2cost_map;

  ErrorMetric(
      HalfedgeGraph &g,
      PointMap &pm,
      KeptPosition< HalfedgeGraph, PointMap > *vkept,
      FEVV::Filters::UniformDequantization< HalfedgeGraph, PointMap > &dequantiz)
      : _g(g), _gt(Geometry(g)), _pm(pm),
        _dequantiz(dequantiz)
  {
    _vkept = vkept;
    _threshold = 0;
  }
  virtual ~ErrorMetric(){}

  virtual void ComputeError() = 0;
  bool isQueueEmpty() { return _queue.empty(); }
  std::tuple< typename boost::graph_traits< HalfedgeGraph >::edge_descriptor,
              double,
              typename Geometry::Point >
  getTopQueue() const { return _queue.top(); }
  void popQueue() { _queue.pop(); }

  double getWeightTop() const { return std::get< 1 >(_queue.top()); }
  virtual double ComputeCostEdge(edge_descriptor e, const Point &collapsePos) = 0;
  double getThreshold() const { return _threshold; }
 
  void delete_from_descriptors(edge_descriptor e) { _edges_cost.erase(e); }

  bool is_present_in_map(edge_descriptor e) const 
  {
    halfedge_descriptor h = halfedge(e, _g);
    halfedge_descriptor h_opp = opposite(h, _g);
    bool is_normal_present = !(_edges_cost.find(edge(h, _g)) == _edges_cost.end());
	if(is_normal_present)
		return true;
    bool is_opposite_present = !(_edges_cost.find(edge(h_opp, _g)) == _edges_cost.end());

    return is_opposite_present; // for most data structure (except LCC) is_normal_present <=> is_opposite_present
  }
  void remove_old_edges(vertex_descriptor vs, vertex_descriptor vt)
  {
    boost::iterator_range<
        CGAL::Halfedge_around_target_iterator< HalfedgeGraph > >
        iterator_range = CGAL::halfedges_around_target(vs, _g);

    for(auto h_around : iterator_range)
    {
      edge_descriptor e = edge(h_around, _g);

      typename std::map< edge_descriptor, std::pair< double, Point > >::iterator
          it = _edges_cost.find(e);
      if(it != _edges_cost.end())
        _edges_cost.erase(it);

      else
      {
        edge_descriptor h_around_2 = edge(opposite(h_around, _g), _g);
        it = _edges_cost.find(h_around_2);
        if(it != _edges_cost.end())
          _edges_cost.erase(it);
        else
          std::cerr << "remove_old_edges: edge not found " << std::endl;
      }
    }


    // for vt
    iterator_range = CGAL::halfedges_around_target(vt, _g);

    for(auto h_around : iterator_range)
    {
      // skip the edge to collapse
      if((source(h_around, _g) == vs && target(h_around, _g) == vt) ||
         (source(h_around, _g) == vt && target(h_around, _g) == vs))
        continue;
      else
      {
        edge_descriptor e = edge(h_around, _g);

        typename std::map< edge_descriptor,
                           std::pair< double, Point > >::iterator it =
            _edges_cost.find(e);
        if(it != _edges_cost.end())
          _edges_cost.erase(it);
        else
        {
          halfedge_descriptor h2 = opposite(h_around, _g);
          edge_descriptor h_around_2 = edge(h2, _g);
          it = _edges_cost.find(h_around_2);
          if(it != _edges_cost.end())
            _edges_cost.erase(it);
          else
            std::cerr << "remove_old_edges: edge not found " << std::endl;
        }
      }
    }
  }

  size_t get_size_queue() const { return _queue.size(); }

  KeptPosition< HalfedgeGraph, PointMap > *
  get_vkept()
  {
    return _vkept;
  }

  virtual std::string getMethodasString() const = 0;

protected:
  HalfedgeGraph &_g;
  const Geometry _gt;
  PointMap &_pm;

  priority_queue_edges _queue; // queue with the weight as the key
  edge2cost_map _edges_cost;

  FEVV::Filters::VKEPT_POSITION _operator;
  KeptPosition< HalfedgeGraph, PointMap >
      *_vkept;
  double _threshold;
  FEVV::Filters::UniformDequantization< HalfedgeGraph, PointMap > &_dequantiz;
};
} // namespace Filters
} // namespace FEVV
