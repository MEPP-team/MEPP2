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
