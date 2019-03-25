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

#include <osg/PolygonMode>
#include <osg/ShapeDrawable>
#include <osg/Material>

#include <osgText/Text>

#include "Base/Color.hpp"

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
           const Color &_color = Color::Clouds(),
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
               const Color &_color = Color::Lime(),
               const std::string &_name = "Cylinder",
               osg::ref_ptr< osg::Group > _group = new osg::Group);

osg::ref_ptr< osg::Group >
createBox(const double &_x,
          const double &_y,
          const double &_z,
          const double &_r,
          const Color &_color = Color::Orange(),
          const std::string &_name = "Box",
          osg::ref_ptr< osg::Group > _group = new osg::Group);

osg::ref_ptr< osg::Group >
createBall(const double &_x,
           const double &_y,
           const double &_z,
           const double &_r,
           const Color &_color = Color::Amethyst(),
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

/**
 * Print OSG tree starting from given node.
 * 
 * \param  nd     node to start from
 * \param  level  first level of indentation, usually 0
 */
void print_osg_tree_from_node(osg::Node *nd, int level = 0);

} // namespace Debug
} // namespace FEVV


// implementation
#include "OSGDebug.inl"
