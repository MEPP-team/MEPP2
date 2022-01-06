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