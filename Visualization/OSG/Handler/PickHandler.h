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

#include <osg/io_utils>
//#include <osgText/Text>
#include "Visualization/SimpleViewer.h"
#include "Visualization/SimpleWindow.h"

//#include <sstream>

namespace FEVV {

// class to handle events with a pick
class PickHandler : public osgGA::GUIEventHandler
{
public:
  PickHandler(FEVV::SimpleViewer *smpViewer,
              osgText::Text *updateText)
      : _smpViewer(smpViewer), _updateText(updateText)
  {
  }

  ~PickHandler() {}

  bool handle(const osgGA::GUIEventAdapter &ea, osgGA::GUIActionAdapter &aa);

  virtual void pick(osgViewer::View *view, const osgGA::GUIEventAdapter &ea);

  void setLabel(const std::string &name)
  {
    if(_updateText.get())
      _updateText->setText(name);
  }

protected:
  osg::ref_ptr< osgText::Text > _updateText;
  FEVV::SimpleViewer *_smpViewer = nullptr;
};


inline
bool
PickHandler::handle(const osgGA::GUIEventAdapter &ea,
                                     osgGA::GUIActionAdapter &aa)
{
  if(ea.getModKeyMask() == osgGA::GUIEventAdapter::MODKEY_SHIFT)
  {
    switch(ea.getEventType())
    {
    case(osgGA::GUIEventAdapter::PUSH):
    {
      osgViewer::View *view = dynamic_cast< osgViewer::View * >(&aa);
      if(view)
        pick(view, ea);
      return false;
    }
    /*case(osgGA::GUIEventAdapter::KEYDOWN):
    {
        if (ea.getKey()=='c')
        {
            osgViewer::View* view = dynamic_cast<osgViewer::View*>(&aa);
            osg::ref_ptr<osgGA::GUIEventAdapter> event = new
    osgGA::GUIEventAdapter(ea); event->setX((ea.getXmin()+ea.getXmax())*0.5);
            event->setY((ea.getYmin()+ea.getYmax())*0.5);
            if (view) pick(view,*event);
        }
        return false;
    }*/
    default:
      return false;
    }
  }

  return false;
}


inline
void
PickHandler::pick(osgViewer::View *view,
                                   const osgGA::GUIEventAdapter &ea)
{
  osgUtil::LineSegmentIntersector::Intersections intersections;

  osg::ref_ptr< osg::Node > node;

  std::string gdlist = "";

  if(view->computeIntersections(ea, intersections))
  {
    for(osgUtil::LineSegmentIntersector::Intersections::iterator hitr =
            intersections.begin();
        hitr != intersections.end();
        ++hitr)
    {
      std::ostringstream os;

      if(!hitr->nodePath.empty() && !(hitr->nodePath.back()->getName().empty()))
      {
        // the geodes are identified by name.
        os << "Object \"" << hitr->nodePath.back()->getName() << "\""
           << std::endl;

        if(hitr->nodePath.back()->getDescription(0) == "MESH")
        {
          std::cout << "Object \"" << hitr->nodePath.back()->getName() << "\""
                    << std::endl;

          node = hitr->drawable->getParent(0);

          // ------------------------------------------------------------------------------------------
          if(!_smpViewer->isNodeSelected(
                 node)) // [PICK] : IMPORTANT test --> see keyword [PICK] in
                        // SimpleWindow.inl
          {
            std::vector< osg::Geode * > geodes =
                _smpViewer->getSelectedGeodes();
            for(unsigned i = 0; i < geodes.size(); i++)
            {
              _smpViewer->setNodeSelected(geodes[i], false);
              // std::cout << "FALSE  : " << geodes[i]->getName() << " - " <<
              // geodes[i] << std::endl;
            }

            _smpViewer->setNodeSelected(node, true);
            // std::cout << "TRUE   : " << node->getName() << " - " << node <<
            // std::endl;

            FEVV::SimpleWindow *sw =
                static_cast< FEVV::SimpleWindow * >(_smpViewer->getWindow());
            sw->update(true);
          }

          // std::cout << std::endl;
          // ------------------------------------------------------------------------------------------

          break; // TEMP : add for only one intersection
        }
      }
      else if(hitr->drawable.valid())
      {
        os << "Object \"" << hitr->drawable->className() << "\"" << std::endl;

        // std::cout<<"Object \""<<hitr->drawable->className()<<"\""<<std::endl;
      }

      os << "        local coords vertex(" << hitr->getLocalIntersectPoint()
         << ")"
         << "  normal(" << hitr->getLocalIntersectNormal() << ")" << std::endl;
      os << "        world coords vertex(" << hitr->getWorldIntersectPoint()
         << ")"
         << "  normal(" << hitr->getWorldIntersectNormal() << ")" << std::endl;

      const osgUtil::LineSegmentIntersector::Intersection::IndexList &vil =
          hitr->indexList;
      for(unsigned int i = 0; i < vil.size(); ++i)
      {
        os << "        vertex indices [" << i << "] = " << vil[i] << std::endl;
      }

      gdlist += os.str();
    }
  }

  setLabel(gdlist);
}

} // namespace FEVV
