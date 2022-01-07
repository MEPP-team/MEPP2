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