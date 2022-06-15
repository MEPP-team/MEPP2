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

#include <QWidget>
#include <QMdiArea>
#include <QTimer>

#if(FEVV_USE_QT5)
#include <osgQOpenGL/osgQOpenGLWidget> // not used ?
#include <osgQOpenGL/osgQOpenGLWindow>
#else
#include <osgQt/GraphicsWindowQt>
// #include <osgViewer/GraphicsWindow>
#endif

#if __GNUC__ >= 9
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#endif

#include <osgViewer/ViewerEventHandlers>

#if __GNUC__ >= 9
#pragma GCC diagnostic pop
#endif

#ifndef Q_MOC_RUN // MT : very important to avoid the error : ' Parse error at
                  // "BOOST_JOIN" ' -> (qt4 pb with boost)
#include "Visualization/BaseAdapterVisuQt.h"
#endif

namespace FEVV {

/**
 * \class SimpleAdapterVisu
 * \brief SimpleAdapterVisu is a specialization of QWidget.
 * This class is the widget added to a QMainWindow, who host an OSG Viewer.
 *
 * @see testViewer.cpp
 */
class SimpleAdapterVisu : public BaseAdapterVisuQt
{
  // Q_OBJECT
public:
  /**
   * Constructor.
   *
   * @note The OpenSceneGraph viewer will be create and initialize if not set
   * (see other constructors).
   *
   * @param[in]   _parent   Pointer to the QWidget parent of this QWidget (used
   * by Qt) (Default value = 0).
   * @param[in]   _f        Windows flags (used by Qt) (Default value = Qt::Widget).
   */
  SimpleAdapterVisu(QWidget *_parent = 0, Qt::WindowFlags _f = Qt::Widget);

  ~SimpleAdapterVisu()
  {
#ifdef DEBUG_VISU2
    std::cout << "*** this=" << this << "    entering " << __func__ << std::endl;
#endif

    // std::cout << "-> ~SimpleAdapterVisu BEGIN" << std::endl;
    delete myViewer;

    // std::cout << "-> ~SimpleAdapterVisu END" << std::endl << std::endl;

#ifdef DEBUG_VISU2
    std::cout << "*** this=" << this << "    leaving " << __func__ << std::endl;
#endif
  }

  void init() override;

  void init(const bool _useMdiWindows);

  bool isInit() const override;

  bool isValid() const override;

  void apply(Plugin *myPlugin) override;

  // void setActive( const bool val, const int freq = 1000 );

  // void attachMesh( HalfedgeGraph& _mesh );

protected:
#if(FEVV_USE_QT5)
  // empty
#else
  virtual void paintEvent(QPaintEvent * /*event*/) override
  {
    static_cast< BaseViewerOSG * >(myViewer)->frame();
  }
#endif

  /**
   * Needed to disable event listener on subwidget.
   * See http://falsinsoft.blogspot.fr/2014/04/qt-get-child-controls-events.html
   * Without this, we can't capture keyPress event.
   **/
  bool eventFilter(QObject *obj, QEvent *event) override;
  //virtual void keyPressEvent(QKeyEvent *event) override; // not used ?
  //virtual void keyReleaseEvent(QKeyEvent *event) override; // not used ?
  virtual bool event(QEvent *event) override;


#if(FEVV_USE_QT5)
  // empty
#else
  void addViewWidget(osg::ref_ptr< osgQt::GraphicsWindowQt > _gw,
                     osg::ref_ptr< osg::Node > _scene);

  osg::ref_ptr< osgQt::GraphicsWindowQt >
  createGraphicsWindow(int _x,
                       int _y,
                       int _w,
                       int _h,
                       const std::string &_name = "",
                       bool _windowDecoration = false) const;
#endif

  // not used ?
  //osgGA::EventQueue *getEventQueue() const;

#if(FEVV_USE_QT5)
public:
  osgQOpenGLWidget *my_osgQOpenGLWidget = nullptr; // not used ?
  osgQOpenGLWindow *my_osgQOpenGLWindow = nullptr;

protected:
#else
protected:
  osgQt::GraphicsWindowQt *myGraphicsWindow = nullptr;
#endif

  QTimer timerUpdate;

  // std::vector< std::pair< HalfedgeGraph, bool > > meshes;
};

} // namespace FEVV

#include "Visualization/SimpleAdapterVisu.inl"
