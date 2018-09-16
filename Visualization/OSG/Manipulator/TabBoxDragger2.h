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
