/**
 * \file		helpers.hxx
 * \author		Vincent Vidal
 * \date     	2017-12-01
 */

#ifndef __EdgeComparators_hxx
#define __EdgeComparators_hxx

#include <limits>

namespace FEVV { 
namespace Comparator 
{
  template <typename Graph, typename PointMap, typename GeometryTraits = FEVV::Geometry_traits<Graph> >
  class EdgeComparator
  {
  public:
    typedef boost::graph_traits<Graph>     GraphTraits;
    typedef typename GraphTraits::vertex_descriptor           vertex_descriptor;
    typedef typename GraphTraits::edge_descriptor             edge_descriptor;
    typedef FEVV::Geometry_traits<Graph>   Geometry;
    typedef typename Geometry::Scalar                         Scalar;
    typedef typename Geometry::Point                          Point;
  private:
    const Graph& _g;
    const PointMap& _pm;
    const GeometryTraits _gt;
  public:
    EdgeComparator(const Graph& g, const PointMap& pm) :_g(g), _pm(pm), _gt(GeometryTraits(g)) {}
    EdgeComparator(const Graph& g, const PointMap& pm, const GeometryTraits& gt):_g(g), _pm(pm), _gt(gt) {}
    bool operator()(edge_descriptor e1, edge_descriptor e2)
    {
      Point pe1_pv1 = get(_pm, source(e1, _g));
      Point pe1_pv2 = get(_pm, target(e1, _g));

      Point pe2_pv1 = get(_pm, source(e2, _g));
      Point pe2_pv2 = get(_pm, target(e2, _g));

      Point pv1((_gt.get_x(pe1_pv1) + _gt.get_x(pe1_pv2))*0.5,
                (_gt.get_y(pe1_pv1) + _gt.get_y(pe1_pv2))*0.5,
                (_gt.get_z(pe1_pv1) + _gt.get_z(pe1_pv2))*0.5);

      Point pv2((_gt.get_x(pe2_pv1) + _gt.get_x(pe2_pv2))*0.5,
                (_gt.get_y(pe2_pv1) + _gt.get_y(pe2_pv2))*0.5,
                (_gt.get_z(pe2_pv1) + _gt.get_z(pe2_pv2))*0.5);

      if (fabs(_gt.get_x(pv1) - _gt.get_x(pv2)) > std::numeric_limits<Scalar>::epsilon())
        return _gt.get_x(pv1) < _gt.get_x(pv2);
      else if (fabs(_gt.get_y(pv1) - _gt.get_y(pv2)) > std::numeric_limits<Scalar>::epsilon())
        return _gt.get_y(pv1) < _gt.get_y(pv2);
      else
        return _gt.get_z(pv1) < _gt.get_z(pv2);
    }
  };
  /////////////////////////////////////////////////////////////////////////////
  template <typename Graph, typename PointMap, typename GeometryTraits = FEVV::Geometry_traits<Graph> >
  static
  EdgeComparator<Graph, PointMap, GeometryTraits>
    get_edge_comparator(const Graph& g, PointMap pm) { return EdgeComparator<Graph, PointMap, GeometryTraits>(g, pm); }
} // namespace Container
} // namespace FEVV

#endif /* __EdgeComparators_hxx */
