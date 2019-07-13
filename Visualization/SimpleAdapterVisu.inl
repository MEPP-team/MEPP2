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
#include "Visualization/BaseViewerOSG.h"

#include "Visualization/Helpers/OSGDebug.hpp"
#include "Visualization/Helpers/QtHelpers.h"

#include <osgGA/MultiTouchTrackballManipulator>

#include <QKeyEvent>
#include <QTextEdit>     // test
#include <QMdiSubWindow> // test

#include <iostream>

#include <functional>
#include <thread>

#include "Visualization/SimpleViewer.h" // @todo temp for getSelectedMeshes()

#include "Visualization/OSG/Handler/PickHandler.h"


inline
FEVV::SimpleAdapterVisu::SimpleAdapterVisu(QWidget *_parent,
                                                            Qt::WindowFlags _f)
    : BaseAdapterVisuQt(_parent, _f)
{
#ifdef DEBUG_VISU2
  std::cout << "*** this=" << this << "    entering " << __func__ << std::endl;
#endif

#ifdef DEBUG_VISU2
  std::cout << "*** this=" << this << "    leaving " << __func__ << std::endl;
#endif
}


inline
void
FEVV::SimpleAdapterVisu::init()
{
  init(false);
}


inline
void
FEVV::SimpleAdapterVisu::init(const bool _useMdiWindows)
{
  if(!Assert::check(
         !bIsInit, "is already init. Leaving...", "SimpleAdapterVisu::init"))
  {
    return;
  }

  if(!Assert::check(myViewer != nullptr,
                    "No viewer attached. Leaving ...",
                    "SimpleAdapterVisu::init"))
  {
    return;
  }

  // installEventFilter( this ); // @FIXME removed to debug child calls.

  connect(&timerUpdate, SIGNAL(timeout()), this, SLOT(update()));
  timerUpdate.start(10);

  // reduce fps when app is idle
  // setActive(false, 100);

#ifdef DEBUG_VISU
  Helpers::changeBackgroundColor(this, Color::Red());
#endif

  // myGraphicsWindow = new osgViewer::GraphicsWindowEmbedded( this->x(),
  // this->y(), this->width(),   this->height() );

  addViewWidget(createGraphicsWindow(0, 0, 500, 500, "Viewer"),
                static_cast< BaseViewerOSG * >(myViewer)->getRootNode());

  // @todo Setting Camera, lights, etc.

  bIsInit = true;
}


inline
bool
FEVV::SimpleAdapterVisu::isInit() const
{
  return bIsInit;
}


inline
bool
FEVV::SimpleAdapterVisu::isValid() const
{
  return bIsInit && (myViewer != nullptr);
}


inline
void
FEVV::SimpleAdapterVisu::addViewWidget(
    osg::ref_ptr< osgQt::GraphicsWindowQt > _gw,
    osg::ref_ptr< osg::Node > _scene)
{
  if(myGraphicsWindow == nullptr)
  {
    myGraphicsWindow = _gw.get();
  }

  // for osgViewer::CompositeViewer
  /*osg::ref_ptr< osgViewer::View > view = new osgViewer::View;
  static_cast< BaseViewerOSG* >( myViewer )->addView( view );*/
  // for osgViewer::Viewer
  BaseViewerOSG *bv = static_cast< BaseViewerOSG * >(myViewer);
  osgViewer::View *view = dynamic_cast< osgViewer::View * >(bv);

  osg::Camera *camera = view->getCamera();
  camera->setGraphicsContext(myGraphicsWindow);

  const osg::GraphicsContext::Traits *traits = myGraphicsWindow->getTraits();

  camera->setClearColor(osg::Vec4(0.2, 0.2, 0.6, 1.0));
  camera->setViewport(new osg::Viewport(0, 0, traits->width, traits->height));

  // set the draw and read buffers up for a double buffered window
  // with rendering going to back buffer
  camera->setDrawBuffer(GL_BACK);
  camera->setReadBuffer(GL_BACK);

  camera->setProjectionMatrixAsPerspective(
      30.0f,
      static_cast< double >(traits->width) /
          static_cast< double >(traits->height),
      1.0f,
      10000.0f);

  // picking
  SimpleViewer *viewer =
      dynamic_cast< SimpleViewer* >(myViewer);
  viewer->hudText = new osgText::Text;
  view->addEventHandler(
      new PickHandler(viewer, viewer->hudText.get()));
  // picking

  view->setSceneData(_scene);
  view->addEventHandler(new osgViewer::StatsHandler);
  view->setCameraManipulator(new osgGA::MultiTouchTrackballManipulator);
  // gw->setTouchEventsEnabled( true );

  auto *osgWidget = myGraphicsWindow->getGLWidget();
  osgWidget->setParent(this);

#ifdef DEBUG_VISU
  Helpers::changeBackgroundColor(osgWidget, Color::Yellow());
#endif

  layout->addWidget(osgWidget);
}


inline
osg::ref_ptr< osgQt::GraphicsWindowQt >
FEVV::SimpleAdapterVisu::createGraphicsWindow(
    int _x,
    int _y,
    int _w,
    int _h,
    const std::string &_name,
    bool _windowDecoration) const
{
  osg::DisplaySettings *ds = osg::DisplaySettings::instance().get();
  osg::ref_ptr< osg::GraphicsContext::Traits > traits =
      new osg::GraphicsContext::Traits;
  traits->windowName = _name;
  traits->windowDecoration = _windowDecoration;
  traits->x = _x;
  traits->y = _y;
  traits->width = _w;
  traits->height = _h;
  traits->doubleBuffer = true;
  traits->alpha = ds->getMinimumNumAlphaBits();
  traits->stencil = ds->getMinimumNumStencilBits();
  traits->sampleBuffers = ds->getMultiSamples();

  //traits->samples = ds->getNumMultiSamples();
  traits->samples = 4; // important : now anti-aliasing is activated here with this value

  traits->vsync = 0; // 0 = disable vsync and 1 = enable vsync

  return new osgQt::GraphicsWindowQt(traits.get());
}

// template< typename HalfedgeGraph >
// void
// FEVV::SimpleAdapterVisu<HalfedgeGraph>::attachMesh( HalfedgeGraph& _mesh )
// {
//     meshes.push_back( {_mesh, false } );
// }

// template< typename HalfedgeGraph >
// typename FEVV::SimpleAdapterVisu<HalfedgeGraph>::Viewer*
// FEVV::SimpleAdapterVisu<HalfedgeGraph>::getViewer()
// {
//     Assert::check( myViewer != nullptr )
//     BOOST_ASSERT_MSG( myViewer != nullptr, "[SimpleAdapterVisu::getViewer()]
//     No viewer attached. Leaving ..." ); if( myViewer == nullptr )
//     {
//         std::cerr << "[SimpleAdapterVisu::getViewer()] No viewer attached.
//         Leaving ..." << std::endl;
//     }
//     return myViewer;
// }

// template< typename HalfedgeGraph >
// typename FEVV::SimpleAdapterVisu<HalfedgeGraph>::Window*
// FEVV::SimpleAdapterVisu<HalfedgeGraph>::getWindow()
// {
//     BOOST_ASSERT_MSG( myWindow != nullptr, "[SimpleAdapterVisu] No window
//     attached. Leaving ..." ); if( myWindow != nullptr )
//     {
//         std::cerr << "[SimpleAdapterVisu] No window attached. Leaving ..." <<
//         std::endl;
//     }
//     return myWindow;
// }

////////////////////////
//// Event handlers ////
////////////////////////


inline
void
FEVV::SimpleAdapterVisu::keyPressEvent(QKeyEvent *event)
{
  QString keyString = event->text();
  const char *keyData = keyString.toLocal8Bit().data();

  if(event->modifiers() == Qt::ControlModifier && event->key() == Qt::Key_S)
  {
    std::cout << "Received CTRL + S " << event->text().toStdString()
              << std::endl;
    myViewer->saveScreenshot("export_screenshot"); // @FIXME Only viewer #0
    return;
  }

  this->getEventQueue()->keyPress(osgGA::GUIEventAdapter::KeySymbol(*keyData));
}


inline
void
FEVV::SimpleAdapterVisu::keyReleaseEvent(QKeyEvent *event)
{
  QString keyString = event->text();
  const char *keyData = keyString.toLocal8Bit().data();

  this->getEventQueue()->keyRelease(
      osgGA::GUIEventAdapter::KeySymbol(*keyData));
}


inline
osgGA::EventQueue *
FEVV::SimpleAdapterVisu::getEventQueue() const
{
  Assert::check(myGraphicsWindow != nullptr,
                "myGraphicsWindow is not initialized.",
                "SimpleAdapterVisu::getEventQueue");

  osgGA::EventQueue *eventQueue = myGraphicsWindow->getEventQueue();

  if(eventQueue)
  {
    return eventQueue;
  }
  else
  {
    throw std::runtime_error("Unable to obtain valid event queue");
  }
}


inline
bool
FEVV::SimpleAdapterVisu::event(QEvent *_event)
{
  bool handled = QWidget::event(_event);

  // if( event->type() != QEvent::Paint)
  // std::cout << "Received event " << event->type() << std::endl;

  // This ensures that the OSG widget is always going to be repainted after the
  // user performed some interaction. Doing this in the event handler ensures
  // that we don't forget about some event and prevents duplicate code.
  switch(_event->type())
  {
  case QEvent::KeyPress:
  case QEvent::KeyRelease:
  case QEvent::ShortcutOverride:
  case QEvent::Resize:
    // case QEvent::MouseButtonDblClick:
    // case QEvent::MouseButtonPress:
    // case QEvent::MouseButtonRelease:
    // case QEvent::MouseMove:
    // case QEvent::Wheel:
    this->update();
    break;

  default:
    break;
  }

  return handled;
}


inline
bool
FEVV::SimpleAdapterVisu::eventFilter(QObject *_obj,
                                                      QEvent *_event)
{
  switch(_event->type())
  {
  // case QEvent::ChildAdded:
  // {
  //     QChildEvent* ce = static_cast<QChildEvent*>( _event );
  //     // Install the filter to each new child object created
  //     // ce->child()->installEventFilter(this); // @FIXME removed to debug
  //     child calls. break;
  // }
  // case QEvent::ChildRemoved:
  // {
  //     QChildEvent* ce = static_cast<QChildEvent*>( _event );
  //     // Remove the the filter from each new child object removed
  //     // ce->child()->removeEventFilter(this); // @FIXME removed to debug
  //     child calls. break;
  // }
  case QEvent::KeyPress:
  case QEvent::KeyRelease:
  {
    QKeyEvent *ke = static_cast< QKeyEvent * >(_event);
    event(ke);

    return true;
  }

  default:
    break;
  }

  return QWidget::eventFilter(_obj, _event);
}

// template< typename HalfedgeGraph >
// std::vector< HalfedgeGraph >
// FEVV::SimpleAdapterVisu<HalfedgeGraph>::getSelectedMeshes()
// {
//     std::vector< HalfedgeGraph> result;
//     for_each( meshes.begin(), meshes.end(),
//     [&result](std::pair<HalfedgeGraph,bool> m)
//     {
//         if(m.second)
//         {
//             result.push_back(m.first);
//         }
//     } );
//     return myViewer->getSelectedMeshes();
// }


inline
void
FEVV::SimpleAdapterVisu::apply(Plugin *myPlugin)
{
  FEVV::MixedMeshesVector meshes;
  std::vector< FEVV::PMapsContainer * > properties_maps;

  auto viewer = dynamic_cast< SimpleViewer * >(myViewer);

  if(viewer)
  {
    meshes = viewer->getSelectedMeshes();
    properties_maps = viewer->getSelected_properties_maps();
  }

  if(meshes.size() > 0)
  {
    // case where a mesh is already loaded and selected
    // apply the plugin to the selected mesh

    for(unsigned i = 0; i < meshes.size(); i++)
    {
      std::cout << "Applying filter" << std::endl;

#ifdef FEVV_USE_CGAL
      if(meshes[i].second == "POLYHEDRON")
      {
        auto mesh_ptr = static_cast< FEVV::MeshPolyhedron* >(meshes[i].first);
        myPlugin->apply(this, mesh_ptr, properties_maps[i]);
      }
      if(meshes[i].second == "SURFACEMESH")
      {
        auto mesh_ptr = static_cast< FEVV::MeshSurface* >(meshes[i].first);
        myPlugin->apply(this, mesh_ptr, properties_maps[i]);
      }
      if(meshes[i].second == "LCC")
      {
        auto mesh_ptr = static_cast< FEVV::MeshLCC* >(meshes[i].first);
        myPlugin->apply(this, mesh_ptr, properties_maps[i]);
      }
      if(meshes[i].second == "CGALPOINTSET")
      {
        auto mesh_ptr = static_cast< FEVV::CGALPointSet* >(meshes[i].first);
        myPlugin->apply(this, mesh_ptr, properties_maps[i]);
      }
#endif //FEVV_USE_CGAL

#ifdef FEVV_USE_OPENMESH
      if(meshes[i].second == "OPENMESH")
      {
        auto mesh_ptr = static_cast< FEVV::MeshOpenMesh* >(meshes[i].first);
        myPlugin->apply(this, mesh_ptr, properties_maps[i]);
      }
#endif //FEVV_USE_OPENMESH

#ifdef FEVV_USE_AIF
      if(meshes[i].second == "AIF")
      {
        auto mesh_ptr = static_cast< FEVV::MeshAIF* >(meshes[i].first);
        myPlugin->apply(this, mesh_ptr, properties_maps[i]);
      }
#endif //FEVV_USE_AIF

#ifdef FEVV_USE_PCL
      if(meshes[i].second == "PCLPOINTCLOUD")
      {
        auto mesh_ptr = static_cast< FEVV::PCLPointCloud* >(meshes[i].first);
        myPlugin->apply(this, mesh_ptr, properties_maps[i]);
      }
#endif //FEVV_USE_PCL
    }
  }
  else
  {
    // case where no viewer is opened or no mesh is loaded
    // apply the plugin to a null mesh
    myPlugin->apply(this, static_cast< void * >(nullptr), nullptr);
  }
}

/*template< typename HalfedgeGraph >
void FEVV::SimpleAdapterVisu<HalfedgeGraph>::setActive( const bool val, const
int freq )
{
    if (val)
    {
        timerUpdate.setInterval(10);
    }
    else
    {
        timerUpdate.setInterval(freq);
    }
}*/
