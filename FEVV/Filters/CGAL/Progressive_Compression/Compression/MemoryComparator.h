#pragma once

#include "FEVV/Filters/CGAL/Progressive_Compression/Compression/CollapseInfo.h"
#include "FEVV/Tools/Comparator/SpanningTreeVertexEdgeComparator.hpp"

namespace FEVV {
namespace Filters {
	
template<typename HalfedgeGraph,
         typename PointMap,
         typename vertex_descriptor = typename boost::graph_traits<HalfedgeGraph>::vertex_descriptor,
         typename Point = typename FEVV::Geometry_traits<HalfedgeGraph>::Point,
         typename Geometry = typename FEVV::Geometry_traits<HalfedgeGraph> >
	class MemoryComparator
	{
	private:
      const FEVV::Comparator::SpanningTreeVertexEdgeComparator<HalfedgeGraph,PointMap> &_st;
	public:
        bool operator()(const FEVV::Filters::CollapseInfo<HalfedgeGraph,PointMap> &mem1,
                        const FEVV::Filters::CollapseInfo<HalfedgeGraph,PointMap> &mem2) const
        {
          return _st.operator()(mem1.get_vkept(), mem2.get_vkept());
        }
        MemoryComparator(
              const FEVV::Comparator::SpanningTreeVertexEdgeComparator< HalfedgeGraph,
                                                                  PointMap >
                  &st)
                      : _st(st)
        {}
        ~MemoryComparator() {}
	};

	template<typename HalfedgeGraph,
             typename PointMap,
             typename vertex_descriptor = typename boost::graph_traits<HalfedgeGraph>::vertex_descriptor,
             typename Geometry = typename FEVV::Geometry_traits<HalfedgeGraph>>
	class VertexSpanComparator
    {
	private:
	  const FEVV::Comparator::SpanningTreeVertexEdgeComparator<HalfedgeGraph,PointMap> &_st;
	public:
	  VertexSpanComparator(
              const FEVV::Comparator::SpanningTreeVertexEdgeComparator< HalfedgeGraph,
                                                                  PointMap >
                  &st)
                      : _st(st)
      {}

      ~VertexSpanComparator() {}

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
    class SimpleVertexComparator
    {
    private:
      const HalfedgeGraph &_g;
      const PointMap &_pm;
      const Geometry _gt;
    public:
      SimpleVertexComparator(const HalfedgeGraph &g, const PointMap &pm)
          : _g(g), _pm(pm), _gt(Geometry(g))
      {}

      ~SimpleVertexComparator() {}

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
}
}