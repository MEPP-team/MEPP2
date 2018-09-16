#include "Visualization/OSG/Manipulator/TabBoxDragger2.h"

#include <osg/ShapeDrawable>
#include <osg/Geometry>
#include <osg/LineWidth>
#include <osg/Quat>

using namespace osgManipulator;

TabBoxDragger2::TabBoxDragger2()
{
  for(int i = 0; i < 6; ++i)
  {
    _planeDraggers.push_back(
        new TabPlaneDragger(1.0f)); // MT : here 20.0f by default
    addChild(_planeDraggers[i].get());
    addDragger(_planeDraggers[i].get());
  }

  {
    _planeDraggers[0]->setMatrix(
        osg::Matrix::translate(osg::Vec3(0.0, 0.5, 0.0)));
  }
  {
    osg::Quat rotation;
    rotation.makeRotate(osg::Vec3(0.0f, -1.0f, 0.0f),
                        osg::Vec3(0.0f, 1.0f, 0.0f));
    _planeDraggers[1]->setMatrix(
        osg::Matrix(rotation) *
        osg::Matrix::translate(osg::Vec3(0.0, -0.5, 0.0)));
  }
  {
    osg::Quat rotation;
    rotation.makeRotate(osg::Vec3(0.0f, 0.0f, 1.0f),
                        osg::Vec3(0.0f, 1.0f, 0.0f));
    _planeDraggers[2]->setMatrix(
        osg::Matrix(rotation) *
        osg::Matrix::translate(osg::Vec3(0.0, 0.0, -0.5)));
  }

  {
    osg::Quat rotation;
    rotation.makeRotate(osg::Vec3(0.0f, 1.0f, 0.0f),
                        osg::Vec3(0.0f, 0.0f, 1.0f));
    _planeDraggers[3]->setMatrix(
        osg::Matrix(rotation) *
        osg::Matrix::translate(osg::Vec3(0.0, 0.0, 0.5)));
  }

  {
    osg::Quat rotation;
    rotation.makeRotate(osg::Vec3(1.0f, 0.0f, 0.0f),
                        osg::Vec3(0.0f, 1.0f, 0.0f));
    _planeDraggers[4]->setMatrix(
        osg::Matrix(rotation) *
        osg::Matrix::translate(osg::Vec3(-0.5, 0.0, 0.0)));
  }

  {
    osg::Quat rotation;
    rotation.makeRotate(osg::Vec3(0.0f, 1.0f, 0.0f),
                        osg::Vec3(1.0f, 0.0f, 0.0f));
    _planeDraggers[5]->setMatrix(
        osg::Matrix(rotation) *
        osg::Matrix::translate(osg::Vec3(0.5, 0.0, 0.0)));
  }

  setParentDragger(getParentDragger());
}

TabBoxDragger2::~TabBoxDragger2() {}

void
TabBoxDragger2::setupDefaultGeometry()
{
  for(unsigned int i = 0; i < _planeDraggers.size(); ++i)
    _planeDraggers[i]->setupDefaultGeometry(false);
}

void
TabBoxDragger2::setPlaneColor(const osg::Vec4 &color)
{
  for(unsigned int i = 0; i < _planeDraggers.size(); ++i)
    _planeDraggers[i]->setPlaneColor(color);
}
