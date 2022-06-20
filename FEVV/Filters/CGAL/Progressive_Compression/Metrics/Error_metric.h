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

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <CGAL/boost/graph/iterator.h>
#include "FEVV/Wrappings/Geometry_traits.h"

#include "FEVV/Filters/CGAL/Progressive_Compression/Operators/Midpoint.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Operators/Halfedge.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Uniform_dequantization.h"

#include <iostream>
#include <vector>
#include <map>
#include <queue>
#include <utility>

namespace FEVV {
namespace Filters {
/// Functor template class used by the priority queue member in
/// the Error_metric class. It permits to store tuples made of
/// edge, cost and point in descending order of their cost. 
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
///        It also manages a priority queue of (edge, cost, collapse position).
template<
    typename HalfedgeGraph,
    typename PointMap>
class Error_metric
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
    priority_queue_edges; /// Priority queue object, manages the priority of
                          /// (edge, cost, kept vertex position) tuples
                          /// according to their cost.

  typedef std::map< edge_descriptor, std::pair< double, Point > > edge2cost_map;

  Error_metric(
      HalfedgeGraph &g,
      PointMap &pm,
      Kept_position< HalfedgeGraph, PointMap > *vkept,
      FEVV::Filters::Uniform_dequantization< HalfedgeGraph, PointMap > &dequantiz)
      : _g(g), _gt(Geometry(g)), _pm(pm),
        _dequantiz(dequantiz)
  {
    _vkept = vkept;
    _threshold = 0;
  }
  virtual ~Error_metric(){}

  /// Method to compute 
  /// 1) all edge costs of the mesh;
  /// 2) the mean cost threshold.
  virtual void compute_error() = 0;
  
  /// Method to compute the cost associated with an edge to collapse.
  /// It usually depends on the resulting vertex position (the collapsed 
  /// position).
  virtual double compute_cost_edge(edge_descriptor e, const Point &collapsePos) = 0;
  
  bool is_queue_empty() { return _queue.empty(); }
  std::tuple< typename boost::graph_traits< HalfedgeGraph >::edge_descriptor,
              double,
              typename Geometry::Point >
  get_top_queue() const { return _queue.top(); }
  void pop_queue() { _queue.pop(); }

  double get_weight_top() const { return std::get< 1 >(_queue.top()); }

  double get_threshold() const { return _threshold; }
 
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

  Kept_position< HalfedgeGraph, PointMap > *
  get_vkept()
  {
    return _vkept;
  }

  virtual std::string get_as_string() const = 0;

protected:
  HalfedgeGraph &_g;
  const Geometry _gt;
  PointMap &_pm;

  priority_queue_edges _queue; /// queue with the weight as the key
  edge2cost_map _edges_cost;

  FEVV::Filters::VKEPT_POSITION _operator;
  Kept_position< HalfedgeGraph, PointMap >
      *_vkept;
  double _threshold;
  FEVV::Filters::Uniform_dequantization< HalfedgeGraph, PointMap > &_dequantiz;
};
} // namespace Filters
} // namespace FEVV
