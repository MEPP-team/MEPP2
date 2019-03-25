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

#include "Visualization/Helpers/OSGHelpers.h"

#include <osg/PolygonMode>
#include <osg/ShapeDrawable>
#include <osg/Material>
#include <osg/Geometry>
#include <osg/Geode>

#include <osg/Texture2D>
#include <osgDB/ReadFile>

inline
osg::ref_ptr< osg::Group >
FEVV::Debug::createLine(const double &_x,
                        const double &_y,
                        const double &_z,
                        const double &_x2,
                        const double &_y2,
                        const double &_z2,
                        const FEVV::Color &_color,
                        const std::string &_name,
                        osg::ref_ptr< osg::Group > _group)
{
  osg::ref_ptr< osg::Vec3Array > points = new osg::Vec3Array;
  osg::ref_ptr< osg::Geometry > geometry = new osg::Geometry;
  osg::ref_ptr< osg::Material > pMaterial = new osg::Material;
  osg::ref_ptr< osg::Geode > geode = new osg::Geode;

  points->push_back(osg::Vec3(_x, _y, _z));
  points->push_back(osg::Vec3(_x2, _y2, _z2));

  geometry->setVertexArray(points.get());
  geometry->addPrimitiveSet(new osg::DrawArrays(GL_LINES, 0, 2));

  geode->addDrawable(geometry);

  pMaterial->setDiffuse(osg::Material::FRONT, Helpers::ColorConverter(_color));
  geode->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
  geode->getOrCreateStateSet()->setAttribute(pMaterial,
                                             osg::StateAttribute::OVERRIDE);

  geode->setName(_name);

  _group->addChild(geode);
  return _group;
}

inline
osg::ref_ptr< osg::Group >
FEVV::Debug::createCylinder(const double &_x,
                            const double &_y,
                            const double &_z,
                            const double &_x2,
                            const double &_y2,
                            const double &_z2,
                            const double &_r,
                            const FEVV::Color &_color,
                            const std::string &_name,
                            osg::ref_ptr< osg::Group > _group)
{
  osg::ref_ptr< osg::Cylinder > cylinder;
  osg::ref_ptr< osg::ShapeDrawable > cylinderDrawable;
  osg::ref_ptr< osg::Material > pMaterial = new osg::Material;
  osg::ref_ptr< osg::Geode > geode = new osg::Geode;

  double height = osg::Vec3(_x - _x2, _y - _y2, _z - _z2).length();
  osg::Vec3 center = osg::Vec3((_x + _x2) / 2, (_y + _y2) / 2, (_z + _z2) / 2);
  osg::Vec3 dir = osg::Vec3(
      0,
      0,
      1); // This is the default direction for the cylinders to face in OpenGL
  osg::Vec3 p = osg::Vec3(_x - _x2, _y - _y2, _z - _z2);
  osg::Vec3 t = dir ^ p;
  double angle = acos((dir * p) / p.length());

  cylinder = new osg::Cylinder(center, _r, height);
  cylinder->setRotation(osg::Quat(angle, osg::Vec3(t.x(), t.y(), t.z())));

  cylinderDrawable = new osg::ShapeDrawable(cylinder);
  geode->addDrawable(cylinderDrawable);

  pMaterial->setDiffuse(osg::Material::FRONT, Helpers::ColorConverter(_color));
  geode->getOrCreateStateSet()->setAttribute(pMaterial,
                                             osg::StateAttribute::OVERRIDE);

  geode->setName(_name);

  _group->addChild(geode);
  return _group;
}

inline
osg::ref_ptr< osg::Group >
FEVV::Debug::createBox(const double &_x,
                       const double &_y,
                       const double &_z,
                       const double &_r,
                       const FEVV::Color &_color,
                       const std::string &_name,
                       osg::ref_ptr< osg::Group > _group)
{
  osg::ref_ptr< osg::Box > box;
  osg::ref_ptr< osg::ShapeDrawable > boxDrawable;
  osg::ref_ptr< osg::Material > pMaterial = new osg::Material;
  osg::ref_ptr< osg::Geode > geode = new osg::Geode;

  box = new osg::Box(osg::Vec3(_x, _y, _z), _r);

  boxDrawable = new osg::ShapeDrawable(box);
  geode->addDrawable(boxDrawable);

  pMaterial->setDiffuse(osg::Material::FRONT, Helpers::ColorConverter(_color));
  geode->getOrCreateStateSet()->setAttribute(pMaterial,
                                             osg::StateAttribute::OVERRIDE);

  geode->setName(_name);

  _group->addChild(geode);
  return _group;
}

inline
osg::ref_ptr< osg::Group >
FEVV::Debug::createBall(const double &_x,
                        const double &_y,
                        const double &_z,
                        const double &_r,
                        const FEVV::Color &_color,
                        const std::string &_name,
                        osg::ref_ptr< osg::Group > _group)
{
  osg::ref_ptr< osg::Sphere > sphere;
  osg::ref_ptr< osg::ShapeDrawable > sphereDrawable;
  osg::ref_ptr< osg::Material > pMaterial = new osg::Material;
  osg::ref_ptr< osg::Geode > geode = new osg::Geode;

  osg::Vec3 center = osg::Vec3(_x, _y, _z);
  sphere = new osg::Sphere(center, _r);

  sphereDrawable = new osg::ShapeDrawable(sphere);
  geode->addDrawable(sphereDrawable);

  pMaterial->setDiffuse(osg::Material::FRONT, Helpers::ColorConverter(_color));
  geode->getOrCreateStateSet()->setAttribute(pMaterial,
                                             osg::StateAttribute::OVERRIDE);

  geode->setName(_name);

  _group->addChild(geode);
  return _group;
}

inline
osg::ref_ptr< osg::Group >
FEVV::Debug::createPyramid(const double &_x,
                           const double &_y,
                           const double &_z,
                           const double &_r,
                           const std::string &_name,
                           osg::ref_ptr< osg::Group > _group)
{
  double r2 = _r / 2.0;

  osg::ref_ptr< osg::Geometry > pyramidGeometry = new osg::Geometry;
  osg::ref_ptr< osg::Geode > geode = new osg::Geode;

  geode->addDrawable(pyramidGeometry);

  osg::ref_ptr< osg::Vec3Array > pyramidVertices = new osg::Vec3Array;
  osg::ref_ptr< osg::Vec3Array > pyramidNormals = new osg::Vec3Array;

  osg::Vec3 p0 = osg::Vec3(_x - r2, _y - r2, _z - r2);
  osg::Vec3 p1 = osg::Vec3(_x + r2, _y - r2, _z - r2);
  osg::Vec3 p2 = osg::Vec3(_x + r2, _y + r2, _z - r2);
  osg::Vec3 p3 = osg::Vec3(_x - r2, _y + r2, _z - r2);
  osg::Vec3 p4 = osg::Vec3(_x, _y, _z + r2);

  pyramidVertices->push_back(p0);
  pyramidVertices->push_back(p1);
  pyramidVertices->push_back(p2);
  pyramidVertices->push_back(p3);
  pyramidVertices->push_back(p4);

  pyramidGeometry->setVertexArray(pyramidVertices);

  osg::ref_ptr< osg::DrawElementsUInt > pyramidBase =
      new osg::DrawElementsUInt(osg::PrimitiveSet::QUADS, 0);
  pyramidBase->push_back(3);
  pyramidBase->push_back(2);
  pyramidBase->push_back(1);
  pyramidBase->push_back(0);
  osg::Vec3 n0 = osg::Vec3(0.0f, 0.0f, -1.0f);
  pyramidNormals->push_back(n0);
  pyramidGeometry->addPrimitiveSet(pyramidBase);

  osg::ref_ptr< osg::DrawElementsUInt > pyramidFaceOne =
      new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES, 0);
  pyramidFaceOne->push_back(0);
  pyramidFaceOne->push_back(1);
  pyramidFaceOne->push_back(4);
  osg::Vec3 n1 = ((p4 - p1) ^ (p0 - p1));
  n1 = n1 / n1.length();
  pyramidNormals->push_back(n1);
  pyramidGeometry->addPrimitiveSet(pyramidFaceOne);

  osg::ref_ptr< osg::DrawElementsUInt > pyramidFaceTwo =
      new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES, 0);
  pyramidFaceTwo->push_back(1);
  pyramidFaceTwo->push_back(2);
  pyramidFaceTwo->push_back(4);
  osg::Vec3 n2 = ((p4 - p2) ^ (p1 - p2));
  n2 = n2 / n2.length();
  pyramidNormals->push_back(n2);
  pyramidGeometry->addPrimitiveSet(pyramidFaceTwo);

  osg::ref_ptr< osg::DrawElementsUInt > pyramidFaceThree =
      new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES, 0);
  pyramidFaceThree->push_back(2);
  pyramidFaceThree->push_back(3);
  pyramidFaceThree->push_back(4);
  osg::Vec3 n3 = ((p4 - p3) ^ (p2 - p3));
  n3 = n3 / n3.length();
  pyramidNormals->push_back(n3);
  pyramidGeometry->addPrimitiveSet(pyramidFaceThree);

  osg::ref_ptr< osg::DrawElementsUInt > pyramidFaceFour =
      new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES, 0);
  pyramidFaceFour->push_back(3);
  pyramidFaceFour->push_back(0);
  pyramidFaceFour->push_back(4);
  osg::Vec3 n4 = ((p4 - p0) ^ (p3 - p0));
  n4 = n4 / n4.length();
  pyramidNormals->push_back(n4);
  pyramidGeometry->addPrimitiveSet(pyramidFaceFour);

  pyramidGeometry->setNormalArray(pyramidNormals,
                                  osg::Array::BIND_PER_PRIMITIVE_SET);

  osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
  colors->push_back(Helpers::ColorConverter(Color::Red()));
  colors->push_back(Helpers::ColorConverter(Color::Green()));
  colors->push_back(Helpers::ColorConverter(Color::Blue()));
  colors->push_back(Helpers::ColorConverter(Color::SunFlower()));
  colors->push_back(Helpers::ColorConverter(Color::Wisteria()));

  pyramidGeometry->setColorArray(colors, osg::Array::BIND_PER_PRIMITIVE_SET);

#if(0) // TEMP
  osg::ref_ptr< osg::Vec2Array > texcoords = new osg::Vec2Array(5);
  (*texcoords)[0].set(0.00f, 0.0f); // tex coord for vertex 0
  (*texcoords)[1].set(0.25f, 0.0f); // tex coord for vertex 1
  (*texcoords)[2].set(0.50f, 0.0f); // ""
  (*texcoords)[3].set(0.75f, 0.0f); // ""
  (*texcoords)[4].set(0.50f, 1.0f); // ""
  pyramidGeometry->setTexCoordArray(0, texcoords);

  osg::ref_ptr< osg::Texture2D > texture = new osg::Texture2D;

  // protect from being optimized away as static state
  texture->setDataVariance(osg::Object::DYNAMIC);

  // load an image by reading a file
  osg::ref_ptr< osg::Image > image =
      osgDB::readImageFile("C:\\tmp\\textures_osg\\texture0.tga");
  if(image)
  {
    // assign the texture to the image we read from file
    texture->setImage(image);

    // ---

    // create a new StateSet with default settings
    osg::ref_ptr< osg::StateSet > stateOne = new osg::StateSet();

    // assign texture unit 0 of our new StateSet to the texture we just created
    // and enable the texture
    stateOne->setTextureAttributeAndModes(0, texture, osg::StateAttribute::ON);

    // associate this state set with the Geode that contains the pyramid
    geode->setStateSet(stateOne);
  }
  else
  {
    std::cout << "Couldn't find texture file." << std::endl;
  }
#endif

  geode->setName(_name);

  _group->addChild(geode);
  return _group;
}

inline
osg::ref_ptr< osg::Group >
FEVV::Debug::createGizmo(osg::ref_ptr< osg::Group > _group)
{
  const int size = 1;
  const double weight = 0.03;

  createCylinder(
      0, 0, 0, size, 0, 0, weight, Color::Red(), "(Gizmo) Cylinder", _group);
  createCylinder(
      0, 0, 0, 0, size, 0, weight, Color::Green(), "(Gizmo) Cylinder", _group);
  createCylinder(
      0, 0, 0, 0, 0, size, weight, Color::Blue(), "(Gizmo) Cylinder", _group);

  if(_group->getName() == "")
  {
    _group->setName("Gizmo");
  }

  return _group;
}

inline
osg::ref_ptr< osg::Group >
FEVV::Debug::createUnitGrid(osg::ref_ptr< osg::Group > _group)
{
  const int size = 10;
  const Color color = Color::Clouds();

  for(int x_axis = -size; x_axis <= size; ++x_axis)
  {
    createLine(
        x_axis, -size, 0, x_axis, size, 0, color, "(UnitGrid) Line", _group);
  }

  for(int y_axis = -size; y_axis <= size; ++y_axis)
  {
    createLine(
        -size, y_axis, 0, size, y_axis, 0, color, "(UnitGrid) Line", _group);
  }

  if(_group->getName() == "")
  {
    _group->setName("UnitGrid");
  }

  return _group;
}

inline
osg::ref_ptr< osg::Group >
FEVV::Debug::createHud(osg::ref_ptr< osgText::Text > updateText,
                       osg::ref_ptr< osg::Group > _group)
{
  std::string timesFont("fonts/times.ttf");

  // turn lighting off for the text and disable depth test to ensure its always
  // ontop.
  osg::Vec3 position(0.0f, 0.0f, 0.0f);

  // this displays what has been selected
  osg::Geode *geode = new osg::Geode();
  osg::StateSet *stateset = geode->getOrCreateStateSet();
  stateset->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
  stateset->setMode(GL_DEPTH_TEST, osg::StateAttribute::OFF);
  geode->setName("(Hud) Hud");
  geode->addDrawable(updateText);

  updateText->setCharacterSize(20.0f);
  updateText->setFont(timesFont);
  updateText->setColor(osg::Vec4(1.0f, 1.0f, 0.0f, 1.0f));
  updateText->setText("");
  updateText->setPosition(position);
  updateText->setDataVariance(osg::Object::DYNAMIC);

  _group->addChild(geode);

  if(_group->getName() == "")
  {
    _group->setName("Hud");
  }

  return _group;
}

inline
void
FEVV::Debug::print_osg_tree_from_node(osg::Node *nd, int level)
{
  std::string blanks(level*2, ' ');

  osg::Geode *geode = dynamic_cast< osg::Geode * >(nd);
  if(geode)
  {
    // the node is a geode
    std::cout << blanks << "geode " << (void *)geode << geode->getName() << std::endl;
  }
  else
  {
    // the node is a group
    osg::Group *gp = dynamic_cast< osg::Group * >(nd);
    if(gp)
    {
      std::cout << blanks << "group " << gp->getName() << std::endl;
      for(unsigned int ic = 0; ic < gp->getNumChildren(); ic++)
        print_osg_tree_from_node(gp->getChild(ic), level + 1);
    }
    else
    {
      std::cout << blanks << "unknown node " << nd << std::endl;
    }
  }
}