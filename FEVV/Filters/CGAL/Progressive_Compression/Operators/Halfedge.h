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

#include "FEVV/Filters/CGAL/Progressive_Compression/Operators/Kept_position.h"

namespace FEVV {
namespace Filters {
template< typename HalfedgeGraph,
          typename PointMap,
          typename Geometry = typename FEVV::Geometry_traits< HalfedgeGraph > >
class Halfedge : public Kept_position< HalfedgeGraph,
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
  Halfedge(HalfedgeGraph &g,
           PointMap &pm)
      : Super_class(g, pm)
  {
    Super_class::_type = VKEPT_POSITION::HALFEDGE;
    Super_class::_reverse = true;
  }
  Halfedge(HalfedgeGraph &g,
           PointMap &pm,
           Geometry &gt)
      : Super_class(g, pm, gt)
  {
    Super_class::_type = VKEPT_POSITION::HALFEDGE;
    Super_class::_reverse = true;
  }
  ~Halfedge(){}
  Point compute_position(
      typename boost::graph_traits< HalfedgeGraph >::edge_descriptor edge) override
  {
    return get(Super_class::_pm, target(edge, Super_class::_g));
  }
  std::string get_as_string() const override { return "halfedge"; }
};

} // namespace Filters
} // namespace FEVV