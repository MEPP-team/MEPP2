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

#include "KeptPosition.h"

namespace FEVV {
namespace Filters {
template< typename HalfedgeGraph,
          typename PointMap,
          typename Geometry = typename FEVV::Geometry_traits< HalfedgeGraph > >
class Midpoint : public KeptPosition< HalfedgeGraph,
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

  Midpoint(HalfedgeGraph &g,
           PointMap &pm)
      : SuperClass(g, pm)
  {
    SuperClass::_type = VKEPT_POSITION::MIDPOINT;
    SuperClass::_reverse = false;
  }
  Midpoint(HalfedgeGraph &g,
           PointMap &pm,
           Geometry &gt)
      : SuperClass(g, pm, gt)
  {
    SuperClass::_type = VKEPT_POSITION::MIDPOINT;
    SuperClass::_reverse = false;
  }

  std::string getMethodasString() const override { return "midpoint"; }
  //compute midpoint of two positions of an edge
  Point ComputePosition(
      typename boost::graph_traits< HalfedgeGraph >::edge_descriptor edge) override
  {
    const Point& source_point = get(SuperClass::_pm, source(edge, SuperClass::_g));
    const Point& target_point = get(SuperClass::_pm, target(edge, SuperClass::_g));
    Point collapsePos = Point((SuperClass::_gt.get_x(source_point) +
                               SuperClass::_gt.get_x(target_point)) /
                                  2,
                              (SuperClass::_gt.get_y(source_point) +
                               SuperClass::_gt.get_y(target_point)) /
                                  2,
                              (SuperClass::_gt.get_z(source_point) +
                               SuperClass::_gt.get_z(target_point)) /
                                  2);


    collapsePos = Point(std::round(SuperClass::_gt.get_x(collapsePos)),
                        std::round(SuperClass::_gt.get_y(collapsePos)),
                        std::round(SuperClass::_gt.get_z(collapsePos)));


    
    return collapsePos;
  }

  ~Midpoint() {}
};
} // namespace Filters
} // namespace FEVV
