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
  * \brief Class used to store information on a collapse operation
  **/
template<
    typename HalfedgeGraph,
    typename PointMap,
    typename vertex_descriptor =
        typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor,
    typename halfedge_descriptor =
        typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor,
    typename face_descriptor =
        typename boost::graph_traits< HalfedgeGraph >::face_descriptor,
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
  vertex_descriptor _vt;
  vertex_descriptor _vs;
  vertex_descriptor
      _v3; /// vertex v3=source(prev(opposite(h))) before the edge collapse
  vertex_descriptor _v4; /// vertex v4=source(prev(h)) before the edge collapse
  Point _pos_v1, _pos_v2;
  std::tuple< bool, bool, bool, bool > _midpoint_rounds;

  Point _pos_v_removed;

  bool _reverse; // reverse the delta or not

  Point _pos_vkept; /// used to store the initial position of vkept in case we
                    /// choose midpoint placement
                    // Vector delta; /// pos(target(h_collapsed)) -
                    // pos(source(h_collapsed))
					
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
  vertex_descriptor get_vt() const { return _vt; }
  vertex_descriptor get_vs() const { return _vs; }
  vertex_descriptor get_v3() const { return _v3; }
  vertex_descriptor get_v4() const { return _v4; }
  const Point& get_pos_v1() const { return _pos_v1; }
  const Point& get_pos_v2() const { return _pos_v2; }
  const Point& get_pos_vkept() const { return _pos_vkept; }
  const std::vector< Vector >& get_error_prediction() const { return _error_prediction; }
  const std::tuple< bool, bool, bool, bool >& get_midpoint_rounds() const { return _midpoint_rounds; }
  
  bool get_reverse() const { return _reverse; }
  
  int get_num_collapse() const { return _num_collapse; }
  void set_num_collapse(int nb) { _num_collapse = nb; }

  void record_error_info(std::tuple< bool, bool, bool, bool > info)
  {
    _midpoint_rounds = info;
  }
  void record_reverse(bool rev) { _reverse = rev; }

  void reverse_bit()
  {
    std::get< 3 >(_midpoint_rounds) = !std::get< 3 >(_midpoint_rounds);
  }
  
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

  void record_v1_v2(const Point& v1, const Point& v2)
  {
    _pos_v1 = v1;
    _pos_v2 = v2;
  }


  void record_pos_removed(halfedge_descriptor h)
  {
    _pos_v_removed = get(_pm, source(h, _g));
  }
  
  void record_vkept(vertex_descriptor vkept)
  {
    _vkept = vkept; // vkept = target(h) = vt
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
  
  void record_error_prediction(const std::vector< Vector >& pred)
  {
    _error_prediction = pred;
  }

  void record_delta(const Point& pos_vkept, VKEPT_POSITION type)
  {
    Geometry gt(_g);
    const Point& pos_vs = _pos_v_removed;
    const Point& pos_vt = pos_vkept;

    if(type == VKEPT_POSITION::MIDPOINT)
    {
      // get delta
      Vector delta(std::move(gt.sub(pos_vt, pos_vs)));

      _error_prediction.push_back(delta);
    }
    else if(type == VKEPT_POSITION::HALFEDGE)
    {
      // get delta
      Vector opp_delta(std::move(gt.sub(pos_vs, pos_vt)));

      _error_prediction.push_back(opp_delta);
    }
    else // full edge collapse
    {
      Vector delta_vt(std::move(gt.sub(pos_vt, pos_vkept)));
      Vector delta_vs(std::move(gt.sub(pos_vs, pos_vkept)));

      _error_prediction.push_back(delta_vt);
      _error_prediction.push_back(delta_vs);
    }
  }

  void record_pos_vkept(vertex_descriptor vkept)
  {
    _pos_vkept = get(_pm, vkept);
  }
};

} // namespace Filters
} // namespace FEVV
