#pragma once

#include "FEVV/Filters/CGAL/Progressive_Compression/Operators/KeptPosition.h"

namespace FEVV {
namespace Filters {
template< typename HalfedgeGraph,
          typename PointMap,
          typename Geometry = typename FEVV::Geometry_traits< HalfedgeGraph > >
class Halfedge : public KeptPosition< HalfedgeGraph,
                                      PointMap,
                                      Geometry >
{
public:
  using Vector = typename Geometry::Vector;
  using Point = typename Geometry::Point;
  typedef KeptPosition< HalfedgeGraph,
                        PointMap,
                        Geometry >
      SuperClass;
  Halfedge(HalfedgeGraph &g,
           PointMap &pm)
      : SuperClass(g, pm)
  {
    SuperClass::_type = VKEPT_POSITION::HALFEDGE;
    SuperClass::_reverse = true;
  }
  Halfedge(HalfedgeGraph &g,
           PointMap &pm,
           Geometry &gt)
      : SuperClass(g, pm, gt)
  {
    SuperClass::_type = VKEPT_POSITION::HALFEDGE;
    SuperClass::_reverse = true;
  }
  ~Halfedge(){}
  Point ComputePosition(
      typename boost::graph_traits< HalfedgeGraph >::edge_descriptor edge) override
  {
    return get(SuperClass::_pm, target(edge, SuperClass::_g));
  }
  std::string getMethodasString() const override { return "halfedge"; }
};

} // namespace Filters
} // namespace FEVV