/**
 * \file		helpers.hxx
 * \author		Vincent Vidal
 * \date     	2017-12-01
 */

#pragma once

#include <limits>

namespace FEVV { 
namespace Comparator 
{
  template <typename Graph, typename PointMap, typename GeometryTraits = FEVV::Geometry_traits<Graph> >
  class VertexComparator
  {
  public: 
    typedef boost::graph_traits<Graph>                        GraphTraits;
    typedef typename GraphTraits::vertex_descriptor           vertex_descriptor;
    typedef FEVV::Geometry_traits<Graph>                      Geometry;
    typedef typename Geometry::Scalar                         Scalar;
    typedef typename Geometry::Point                          Point;
  private:
    const PointMap& _pm;
    const GeometryTraits _gt;
  public:
    VertexComparator(const Graph& g, const PointMap& pm) :_pm(pm), _gt(GeometryTraits(g)) {}
    VertexComparator(const PointMap& pm, const GeometryTraits& gt) :_pm(pm), _gt(gt) {}
    bool operator()(vertex_descriptor v1, vertex_descriptor v2)
    {
      Point pv1 = get(_pm, v1);
      Point pv2 = get(_pm, v2);

      if (fabs(_gt.get_x(pv1) - _gt.get_x(pv2)) > std::numeric_limits<Scalar>::epsilon())
        return _gt.get_x(pv1) < _gt.get_x(pv2);
      else if(fabs(_gt.get_y(pv1) - _gt.get_y(pv2)) > std::numeric_limits<Scalar>::epsilon())
        return _gt.get_y(pv1) < _gt.get_y(pv2);
      else
        return _gt.get_z(pv1) < _gt.get_z(pv2);
    }
  };
  /////////////////////////////////////////////////////////////////////////////
  template <typename Graph, typename PointMap, typename GeometryTraits = FEVV::Geometry_traits<Graph> >
  static 
  VertexComparator<Graph, PointMap, GeometryTraits>
    get_vertex_comparator(const Graph& g, PointMap pm) { return VertexComparator<Graph, PointMap, GeometryTraits>(g, pm); }
} // namespace Container
} // namespace FEVV

