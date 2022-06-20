// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
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

#include "FEVV/Wrappings/Geometry_traits.h"

#include <limits>
#include <iterator>

namespace FEVV { 
namespace Comparator 
{
  template <typename Graph, typename PointMap, typename GeometryTraits = FEVV::Geometry_traits<Graph> >
  class Vertex_comparator
  {
  public: 
    typedef boost::graph_traits<Graph>                        GraphTraits;
    typedef typename GraphTraits::vertex_descriptor           vertex_descriptor;
    typedef typename GeometryTraits::Scalar                   Scalar;
    typedef typename GeometryTraits::Point                    Point;
  private:
    const Graph& _g;
    const PointMap& _pm;
    const GeometryTraits _gt;
  public:
    Vertex_comparator(const Graph& g, const PointMap& pm) :_g(g), _pm(pm), _gt(GeometryTraits(g)) {}
    Vertex_comparator(const Graph& g, const PointMap& pm, const GeometryTraits& gt) :_g(g), _pm(pm), _gt(gt) {}
    bool operator()(vertex_descriptor v1, vertex_descriptor v2)
    {
      if(v1==v2)
        return false;

      if(v1==GraphTraits::null_vertex())
      {
		  return false;
	  }
	  
      if(v2==GraphTraits::null_vertex())
      {
		  return false;
	  }	  
	  
      Point pv1 = get(_pm, v1);
      Point pv2 = get(_pm, v2);

      if (fabs(_gt.get_x(pv1) - _gt.get_x(pv2)) > std::numeric_limits<Scalar>::epsilon())
        return _gt.get_x(pv1) < _gt.get_x(pv2);
      else if(fabs(_gt.get_y(pv1) - _gt.get_y(pv2)) > std::numeric_limits<Scalar>::epsilon())
        return _gt.get_y(pv1) < _gt.get_y(pv2);
      else if(fabs(_gt.get_z(pv1) - _gt.get_z(pv2)) > std::numeric_limits<Scalar>::epsilon())
        return _gt.get_z(pv1) < _gt.get_z(pv2);
	
      auto e_v1 = in_edges(v1,_g), 
           e_v2 = in_edges(v2,_g);
	     
      auto e_v1_deg = std::distance(e_v1.first, e_v1.second),
           e_v2_deg = std::distance(e_v2.first, e_v2.second);
      if( e_v1_deg != e_v2_deg)
		  return e_v1_deg < e_v2_deg ;
	  
      e_v1_deg=e_v2_deg=0;
      decltype(e_v1_deg) e_v1_min_deg = 1000, e_v1_max_deg = 0, e_v2_min_deg = 1000, e_v2_max_deg = 0;
      auto iter_e_v = e_v1.first;
	  for( ; iter_e_v!=e_v1.second; ++iter_e_v)
	  {
        auto e_vopp = (source(*iter_e_v, _g)!=v1)?in_edges(source(*iter_e_v, _g),_g):in_edges(target(*iter_e_v, _g),_g);
        auto current_deg = std::distance(e_vopp.first, e_vopp.second);
        if (current_deg < e_v1_min_deg)
          e_v1_min_deg = current_deg;
        if (current_deg > e_v1_max_deg)
          e_v1_max_deg = current_deg;
        e_v1_deg += current_deg;
	  }
	  
	  iter_e_v = e_v2.first;
	  for( ; iter_e_v!=e_v2.second; ++iter_e_v)
	  {
        auto e_vopp = (source(*iter_e_v, _g)!=v2)?in_edges(source(*iter_e_v, _g),_g):in_edges(target(*iter_e_v, _g),_g);
        auto current_deg = std::distance(e_vopp.first, e_vopp.second);
        if (current_deg < e_v2_min_deg)
          e_v2_min_deg = current_deg;
        if (current_deg > e_v2_max_deg)
          e_v2_max_deg = current_deg;
        e_v2_deg += current_deg;
	  }
	  
	  if( e_v1_deg != e_v2_deg)
        return e_v1_deg < e_v2_deg ;
	  
      if (e_v1_min_deg != e_v2_min_deg)
        return e_v1_min_deg < e_v2_min_deg;

      if (e_v1_max_deg != e_v2_max_deg)
        return e_v1_max_deg < e_v2_max_deg;

	  return false;
    }
  };
  /////////////////////////////////////////////////////////////////////////////
  template <typename Graph, typename PointMap, typename GeometryTraits = FEVV::Geometry_traits<Graph> >
  static 
  Vertex_comparator<Graph, PointMap, GeometryTraits>
    get_vertex_comparator(const Graph& g, const PointMap& pm) { return Vertex_comparator<Graph, PointMap, GeometryTraits>(g, pm); }
} // namespace Container
} // namespace FEVV

