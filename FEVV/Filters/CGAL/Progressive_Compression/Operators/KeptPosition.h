#pragma once

#include "FEVV/Filters/CGAL/Progressive_Compression/Parameters.h"

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <CGAL/boost/graph/iterator.h>
#include "FEVV/Wrappings/Geometry_traits.h"

namespace FEVV {
namespace Filters {
template< typename HalfedgeGraph,
          typename PointMap,
          typename Geometry = typename FEVV::Geometry_traits< HalfedgeGraph > >
class KeptPosition
{
public:
  using Vector = typename Geometry::Vector;
  using Point = typename Geometry::Point;
  KeptPosition(HalfedgeGraph &g,
               PointMap &pm)
      : _g(g), _gt(Geometry(_g)), _pm(pm){}
  KeptPosition(HalfedgeGraph &g,
               PointMap &pm,
               Geometry &gt)
      : _g(g), _gt(gt), _pm(pm){}
  virtual ~KeptPosition(){}
  virtual Point ComputePosition(
      typename boost::graph_traits< HalfedgeGraph >::edge_descriptor edge) = 0;
  VKEPT_POSITION getType() const { return _type; }
  bool getReverse() const { return _reverse; }
  virtual std::string getMethodasString() const = 0;

protected:
  HalfedgeGraph &_g;
  const Geometry _gt;
  PointMap &_pm;
  FEVV::Filters::VKEPT_POSITION _type;
  bool _reverse;
};
} // namespace Filters
} // namespace FEVV