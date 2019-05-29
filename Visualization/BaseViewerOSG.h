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

#include "Visualization/BaseViewer.h"

#include <osg/Version>

//#include <osgViewer/CompositeViewer>
#include <osgViewer/Viewer>

#include <osgText/Text>

#include "Base/Color.hpp"
#include "Visualization/OSG/Visitor/DataVisitor.h"

// Generic mesh iterators
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include "FEVV/Wrappings/Geometry_traits.h"

namespace FEVV {

class BaseViewerOSG : public osgViewer::Viewer, public BaseViewer
{
public:
  using DataModel = DataVisitor::Data;
  using DataModelVector = DataVisitor::Output;
  using Model = osg::Geode;
  using Group = osg::Group;

  using BaseViewer::isSelected;
  using BaseViewer::setSelected;

public:
public:
  /**
   * Constructor.
   */
  BaseViewerOSG()
      : root_node(new osg::Group),
        visitor(new DataVisitor(this)), osgViewer::Viewer(), BaseViewer()
  {
#ifdef DEBUG_VISU2
    std::cout << "*** this=" << this << "    entering " << __func__ << std::endl;
#endif

#if(FEVV_USE_QT5)
    // Qt5 is currently crashing and reporting "Cannot make QOpenGLContext
    // current in a different thread" when the viewer is run multi-threaded,
    // this is regression from Qt4
#if OSG_MIN_VERSION_REQUIRED(3, 4, 0)
    // setThreadingModel(osgViewer::ViewerBase::AutomaticSelection); // maybe
    // one day this PB will be fixed either by OSG or by Qt...
    setThreadingModel(osgViewer::ViewerBase::SingleThreaded);
#else
    setThreadingModel(osgViewer::ViewerBase::SingleThreaded);
#endif
#else
#if defined(Q_WS_X11)    // X Window System -> Linux
    // setThreadingModel(osgViewer::ViewerBase::CullThreadPerCameraDrawThreadPerContext);
    // // PB with osgFX/Outline -> video driver freeze/crash
    // setThreadingModel(osgViewer::ViewerBase::AutomaticSelection); // PB with
    // some intel graphic cards
    setThreadingModel(osgViewer::ViewerBase::SingleThreaded);
#elif defined(Q_WS_MACX) // Mac OS X
    // setThreadingModel(osgViewer::ViewerBase::CullThreadPerCameraDrawThreadPerContext);
    // // PB with osgFX/Outline -> video driver freeze/crash
    setThreadingModel(osgViewer::ViewerBase::AutomaticSelection);
#else                    // Windows
    // setThreadingModel(osgViewer::ViewerBase::CullDrawThreadPerContext);
    setThreadingModel(osgViewer::ViewerBase::AutomaticSelection);
#endif
#endif

    root_node->setName("Root");

#ifdef DEBUG_VISU2
    std::cout << "*** this=" << this << "    leaving " << __func__ << std::endl;
#endif
  }

  virtual ~BaseViewerOSG()
  {
#ifdef DEBUG_VISU2
    std::cout << "*** this=" << this << "    entering " << __func__ << std::endl;
#endif

    delete visitor;
    root_node->removeChildren(0, root_node->getNumChildren());
    //delete root_node;
    //ELO-note: osg::Group destructor is protected ;
    //          can not be called directly!
    //          error: calling a protected destructor of class 'osg::Group'

#ifdef DEBUG_VISU2
    std::cout << "*** this=" << this << "    leaving " << __func__ << std::endl;
#endif
  }

  /**
   * Get the root node of the scene.
   *
   * @return the root node of the scene.
   */
  Group *getRootNode() { return root_node; }

  /**
   * Add a geode to the scene.
   *
   * @note A geode is a Geometry Node.
   *
   * @param[in] _geode Pointer to a geode.
   **/
  virtual void addModel(Model *_geode) = 0;

  /**
   * Add a group to the scene.
   *
   * @note A group is a set of geodes.
   *
   * @param[in] _group Pointer to a group of geode.
   **/
  virtual void addGroup(Group *_group) = 0;

  virtual void setNodeSelected(osg::Node *_node, bool isSelected) = 0;
  virtual bool isNodeSelected(osg::Node *_node) = 0;

  virtual DataModelVector *getDataModel() = 0;

protected:
  osg::Group *root_node = nullptr;
  DataVisitor *visitor = nullptr;
};

} // namespace FEVV
