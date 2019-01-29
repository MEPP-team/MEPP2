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

#include <osgManipulator/TabPlaneDragger>

namespace osgManipulator {

/**
 * TabBoxDragger2 consists of 6 TabPlaneDraggers to form a box dragger that
 * performs translation and scaling.
 */
class TabBoxDragger2 : public CompositeDragger
{
public:
  TabBoxDragger2();

  META_OSGMANIPULATOR_Object(osgManipulator, TabBoxDragger2)

      /** Setup default geometry for dragger. */
      void setupDefaultGeometry();

  void setPlaneColor(const osg::Vec4 &color);

protected:
  virtual ~TabBoxDragger2();

  std::vector< osg::ref_ptr< TabPlaneDragger > > _planeDraggers;
};

} // namespace osgManipulator
