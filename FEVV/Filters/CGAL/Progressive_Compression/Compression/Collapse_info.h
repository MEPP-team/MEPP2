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
#include "FEVV/Wrappings/Geometry_traits.h"
#include <boost/graph/properties.hpp>
#include "FEVV/Filters/CGAL/Progressive_Compression/Parameters.h"
#include <vector>
#include <map>

namespace FEVV {
namespace Filters {
/** 
  * \brief Class used to store information on a single edge
  *        collapse operation.
  **/
template<
    typename HalfedgeGraph,
    typename PointMap,
    typename vertex_descriptor =
        typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor,
    typename halfedge_descriptor =
        typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor,
    typename Vector = typename FEVV::Geometry_traits< HalfedgeGraph >::Vector,
    typename Point = typename FEVV::Geometry_traits< HalfedgeGraph >::Point,
    typename Geometry = typename FEVV::Geometry_traits< HalfedgeGraph > >
class Collapse_info
{
private:
  HalfedgeGraph &_g;
  const Geometry _gt;
  PointMap &_pm;
  vertex_descriptor _vkept; /// kept vertex after applying the edge collapse
                            /// (the remaining vertex after the operation, 
                            /// i.e. either source or target vertex)
  vertex_descriptor _vt; /// target vertex
  vertex_descriptor _vs; /// source vertex
  vertex_descriptor
      _v3; /// vertex v3=source(prev(opposite(h))) before the edge collapse 
           /// (pivot vertex when not null)
  vertex_descriptor _v4; /// vertex v4=source(prev(h)) before the edge collapse
                         /// (pivot vertex when not null)  
  Point _pos_vt, _pos_vs;

  bool _reverse; /// reverse the delta or not

  Point _pos_vkept; /// used to store the initial position of vkept in case we
                    /// choose midpoint placement
				
  std::vector< Vector >
      _error_prediction; /// Contains one or several deltas for position					
  /////////////////////////////////////////////////////////////////////////////
  int _num_collapse;
public:
  Collapse_info(HalfedgeGraph &g, PointMap &pm) : _g(g), _gt(Geometry(_g)), _pm(pm)
  {
    _reverse = false;
  }
  ~Collapse_info() {}

  vertex_descriptor get_vkept() const { return _vkept; }
  
  /// Get edge target vertex (sometimes called v1).
  vertex_descriptor get_vt() const { return _vt; }
  
  /// Get edge source vertex (sometimes called v0 or v2).
  vertex_descriptor get_vs() const { return _vs; }
  
  /// Get the vertex in front of the edge to collapse and 
  /// on the other triangle.
  /// When not a null_vertex, this vertex is a pivot.
  vertex_descriptor get_v3() const { return _v3; }
  
  /// Get the vertex in front of the edge to collapse and  
  /// on the same triangle.
  /// When not a null_vertex, this vertex is a pivot.  
  vertex_descriptor get_v4() const { return _v4; }
  
  const Point& get_pos_vt() const { return _pos_vt; }
  const Point& get_pos_vs() const { return _pos_vs; }
  const Point& get_pos_vkept() const { return _pos_vkept; }
  const std::vector< Vector >& get_error_prediction() const { return _error_prediction; }
  
  /// Get the reverse information (used by the Predictor object,
  /// see compute_residuals method)
  bool get_reverse() const { return _reverse; }
  
  int get_num_collapse() const { return _num_collapse; }
  void set_num_collapse(int nb) { _num_collapse = nb; }

  /// Set reverse information to true when the vertex v3 is reached
  /// after v4 when going throught adjacent vertices to the kept
  /// vertex. See connectivity_encoding method for more details. 
  void record_reverse(bool rev) { _reverse = rev; }

  void record_v3_v4(halfedge_descriptor h)
  {
    _v3 = boost::graph_traits< HalfedgeGraph >::null_vertex();
    _v4 = boost::graph_traits< HalfedgeGraph >::null_vertex();
    if(!CGAL::is_border(h, _g))
    {
      _v4 = source(prev(h, _g), _g);
	}

    if(!CGAL::is_border(opposite(h, _g), _g))
    {
      _v3 = source(prev(opposite(h, _g), _g), _g);
    }
  }

  void record_vt_vs(const Point& v1, const Point& v2)
  {
    _pos_vt = v1;
    _pos_vs = v2;
  }

  void record_vkept(vertex_descriptor vkept)
  {
    _vkept = vkept; 
  }

  bool is_vkept_source()
  {
    if(_vkept == _vs)
      return true;
    else
      return false;
  }

  void set_source_and_target(vertex_descriptor s, vertex_descriptor t)
  {
    _vs = s;
    _vt = t;
  }
  
  /// Store the computed residuals which is done during the prediction
  /// (see compute_residuals method of any Predictor object).  
  void record_error_prediction(const std::vector< Vector >& pred)
  {
    _error_prediction = pred;
  }

  void record_pos_vkept(vertex_descriptor vkept)
  {
    _pos_vkept = get(_pm, vkept);
  }
};

} // namespace Filters
} // namespace FEVV
