#pragma once

#if defined(OSGDebug_RECURSES)
#error Recursive header files inclusion detected in OSGDebug.h
#else // defined(OSGDebug_RECURSES)
/** Prevents recursive inclusion of headers. */
#define OSGDebug_RECURSES

#if !defined OSGDebug_h
/** Prevents repeated inclusion of headers. */
#define OSGDebug_h

#include <osg/PolygonMode>
#include <osg/ShapeDrawable>
#include <osg/Material>

#include <osgText/Text>

#include "Base/Color.h"

// #define DEBUG_VISU //@FIXME to be removed. Only for GUI debug purpose.


namespace FEVV {

namespace Debug {

osg::ref_ptr< osg::Group >
createLine(const double &_x,
           const double &_y,
           const double &_z,
           const double &_x2,
           const double &_y2,
           const double &_z2,
           const Color &_color = Color::Clouds,
           const std::string &_name = "Line",
           osg::ref_ptr< osg::Group > _group = new osg::Group);

osg::ref_ptr< osg::Group >
createCylinder(const double &_x,
               const double &_y,
               const double &_z,
               const double &_x2,
               const double &_y2,
               const double &_z2,
               const double &_r,
               const Color &_color = Color::Lime,
               const std::string &_name = "Cylinder",
               osg::ref_ptr< osg::Group > _group = new osg::Group);

osg::ref_ptr< osg::Group >
createBox(const double &_x,
          const double &_y,
          const double &_z,
          const double &_r,
          const Color &_color = Color::Orange,
          const std::string &_name = "Box",
          osg::ref_ptr< osg::Group > _group = new osg::Group);

osg::ref_ptr< osg::Group >
createBall(const double &_x,
           const double &_y,
           const double &_z,
           const double &_r,
           const Color &_color = Color::Amethyst,
           const std::string &_name = "Ball",
           osg::ref_ptr< osg::Group > _group = new osg::Group);

osg::ref_ptr< osg::Group >
createPyramid(const double &_x,
              const double &_y,
              const double &_z,
              const double &_r,
              const std::string &_name = "Pyramid",
              osg::ref_ptr< osg::Group > _group = new osg::Group);

osg::ref_ptr< osg::Group >
createGizmo(osg::ref_ptr< osg::Group > _group = new osg::Group);

osg::ref_ptr< osg::Group >
createUnitGrid(osg::ref_ptr< osg::Group > _group = new osg::Group);

osg::ref_ptr< osg::Group >
createHud(osg::ref_ptr< osgText::Text > updateText,
          osg::ref_ptr< osg::Group > _group = new osg::Group);

} // namespace Debug

} // namespace FEVV

#endif // !defined OSGDebug_h

#undef OSGDebug_RECURSES
#endif // else defined(OSGDebug_RECURSES)
