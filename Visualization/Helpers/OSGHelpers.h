// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published 
// by the Free Software Foundation; either version 3 of the License, 
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#pragma once

#include <osg/Vec3>
#include <osg/Vec4>
#include <osg/Group>
#include <osg/Geode>
#include <osg/Geometry>
#include <osg/Material>

#include "Base/Color.hpp"
#include "Visualization/BaseViewer.h"

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include "FEVV/Wrappings/Geometry_traits.h"

namespace FEVV {

namespace Helpers {

static unsigned int nbMeshDrawed = 0;

osg::Vec4
ColorConverter(const Color &_color);

template<
    typename HalfedgeGraph,
    typename Vector = typename FEVV::Geometry_traits< HalfedgeGraph >::Vector >
osg::Vec4
VectorColorConverter(const Vector &_color)
{
  osg::Vec4 color((float)_color[0], (float)_color[1], (float)_color[2], 1.0f);

  return color;
}

template<
    typename HalfedgeGraph,
    typename Vector = typename FEVV::Geometry_traits< HalfedgeGraph >::Vector >
osg::Vec3
VectorConverter(const Vector &_vector)
{
  osg::Vec3 vec(_vector[0], _vector[1], _vector[2]);

  return vec;
}

inline osg::Vec4
ColorConverter(const Color &_color)
{
  osg::Vec4 color((float)_color.red() / 255.0f,
                  (float)_color.green() / 255.0f,
                  (float)_color.blue() / 255.0f,
                  (float)_color.alpha() / 255.0f);

  return color;
}

} // namespace Helpers

} // namespace FEVV
