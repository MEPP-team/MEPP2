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

#include "FEVV/Filters/CGAL/Progressive_Compression/Compression/Collapse_info.h"
#include "FEVV/Tools/Comparator/SpanningTreeVertexEdgeComparator.hpp"

namespace FEVV {
namespace Filters {
	
template<typename HalfedgeGraph,
         typename PointMap,
         typename vertex_descriptor = typename boost::graph_traits<HalfedgeGraph>::vertex_descriptor,
         typename Point = typename FEVV::Geometry_traits<HalfedgeGraph>::Point,
         typename Geometry = typename FEVV::Geometry_traits<HalfedgeGraph> >
class Memory_comparator
{
	private:
      const FEVV::Comparator::SpanningTreeVertexEdgeComparator<HalfedgeGraph,PointMap> &_st;
	public:
        bool operator()(const FEVV::Filters::Collapse_info<HalfedgeGraph,PointMap> &mem1,
                        const FEVV::Filters::Collapse_info<HalfedgeGraph,PointMap> &mem2) const
        {
          return _st.operator()(mem1.get_vkept(), mem2.get_vkept());
        }
        Memory_comparator(
              const FEVV::Comparator::SpanningTreeVertexEdgeComparator< HalfedgeGraph,
                                                                  PointMap >
                  &st)
                      : _st(st)
        {}
        ~Memory_comparator() {}
};

template<typename HalfedgeGraph,
         typename PointMap,
         typename vertex_descriptor = typename boost::graph_traits<HalfedgeGraph>::vertex_descriptor,
         typename Geometry = typename FEVV::Geometry_traits<HalfedgeGraph>>
class Vertex_span_comparator
{
	private:
	  const FEVV::Comparator::SpanningTreeVertexEdgeComparator<HalfedgeGraph,PointMap> &_st;
	public:
	  Vertex_span_comparator(
              const FEVV::Comparator::SpanningTreeVertexEdgeComparator< HalfedgeGraph,
                                                                  PointMap >
                  &st)
                      : _st(st)
      {}

      ~Vertex_span_comparator() {}

      bool operator()(vertex_descriptor v1, vertex_descriptor v2) const
      {
        return _st.operator()(v1, v2);
      }
	
};

template<
        typename HalfedgeGraph,
        typename PointMap,
        typename vertex_descriptor =
            typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor,
        typename Point = typename FEVV::Geometry_traits< HalfedgeGraph >::Point,
        typename Geometry = typename FEVV::Geometry_traits< HalfedgeGraph > >
class Vertex_lexicographic_comparator
{
    private:
      const HalfedgeGraph &_g;
      const PointMap &_pm;
      const Geometry _gt;
    public:
      Vertex_lexicographic_comparator(const HalfedgeGraph &g, const PointMap &pm)
          : _g(g), _pm(pm), _gt(Geometry(g))
      {}

      ~Vertex_lexicographic_comparator() {}

      bool operator()(vertex_descriptor v1, vertex_descriptor v2) const
      {
        const Point& pv1 = get(_pm, v1);
        const Point& pv2 = get(_pm, v2);

        if(fabs(_gt.get_x(pv1) - _gt.get_x(pv2)) >
           std::numeric_limits< double >::epsilon())
          return _gt.get_x(pv1) < _gt.get_x(pv2);
        else if(fabs(_gt.get_y(pv1) - _gt.get_y(pv2)) >
                std::numeric_limits<double >::epsilon())
          return _gt.get_y(pv1) < _gt.get_y(pv2);
        else if(fabs(_gt.get_z(pv1) - _gt.get_z(pv2)) >
                std::numeric_limits< double >::epsilon())
          return _gt.get_z(pv1) < _gt.get_z(pv2);
		  
		return false;
      }
};
} // namespace Filters
} // namespace FEVV