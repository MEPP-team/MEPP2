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

#include "Kept_position.h"

namespace FEVV {
namespace Filters {
/**
 *  \brief Concrete class to represent the midpoint position type of
 *         the resulting vertex of an edge collapse. 
 */	
template< typename HalfedgeGraph,
          typename PointMap,
          typename Geometry = typename FEVV::Geometry_traits< HalfedgeGraph > >
class Midpoint : public Kept_position< HalfedgeGraph,
                                      PointMap,
                                      Geometry >
{
public:
  using Vector = typename Geometry::Vector;
  using Point = typename Geometry::Point;
  typedef Kept_position< HalfedgeGraph,
                        PointMap,
                        Geometry >
      Super_class;

  Midpoint(HalfedgeGraph &g,
           PointMap &pm)
      : Super_class(g, pm)
  {
    Super_class::_type = VKEPT_POSITION::MIDPOINT;
    Super_class::_reverse = false; // has no reverse case
  }
  Midpoint(HalfedgeGraph &g,
           PointMap &pm,
           Geometry &gt)
      : Super_class(g, pm, gt)
  {
    Super_class::_type = VKEPT_POSITION::MIDPOINT;
    Super_class::_reverse = false; // has no reverse case
  }

  std::string get_as_string() const override { return "midpoint"; }
  
  /// Compute the midpoint of an edge.
  Point compute_position(
      typename boost::graph_traits< HalfedgeGraph >::edge_descriptor edge) override
  {
    const Point& source_point = get(Super_class::_pm, source(edge, Super_class::_g));
    const Point& target_point = get(Super_class::_pm, target(edge, Super_class::_g));
    Point collapsePos = Point((Super_class::_gt.get_x(source_point) +
                               Super_class::_gt.get_x(target_point)) /
                                  2,
                              (Super_class::_gt.get_y(source_point) +
                               Super_class::_gt.get_y(target_point)) /
                                  2,
                              (Super_class::_gt.get_z(source_point) +
                               Super_class::_gt.get_z(target_point)) /
                                  2);


    collapsePos = Point(std::round(Super_class::_gt.get_x(collapsePos)),
                        std::round(Super_class::_gt.get_y(collapsePos)),
                        std::round(Super_class::_gt.get_z(collapsePos)));


    
    return collapsePos;
  }

  ~Midpoint() {}
};
} // namespace Filters
} // namespace FEVV
