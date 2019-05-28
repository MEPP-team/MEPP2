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

#include <osgQt/GraphicsWindowQt>
// #include <osgViewer/GraphicsWindow>
#include <osgViewer/ViewerEventHandlers>

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
   * @param[in]   _f        Windows flags (used by Qt) (Default value = 0).
   */
  SimpleAdapterVisu(QWidget *_parent = 0, Qt::WindowFlags _f = 0);

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
  virtual void paintEvent(QPaintEvent *event) override
  {
    static_cast< BaseViewerOSG * >(myViewer)->frame();
  }

  /**
   * Needed to disable event listener on subwidget.
   * See http://falsinsoft.blogspot.fr/2014/04/qt-get-child-controls-events.html
   * Without this, we can't capture keyPress event.
   **/
  bool eventFilter(QObject *obj, QEvent *event) override;
  virtual void keyPressEvent(QKeyEvent *event) override;
  virtual void keyReleaseEvent(QKeyEvent *event) override;
  virtual bool event(QEvent *event) override;

  void addViewWidget(osg::ref_ptr< osgQt::GraphicsWindowQt > _gw,
                     osg::ref_ptr< osg::Node > _scene);
  osg::ref_ptr< osgQt::GraphicsWindowQt >
  createGraphicsWindow(int _x,
                       int _y,
                       int _w,
                       int _h,
                       const std::string &_name = "",
                       bool _windowDecoration = false) const;
  osgGA::EventQueue *getEventQueue() const;

protected:
  osgQt::GraphicsWindowQt *myGraphicsWindow = nullptr;

  QTimer timerUpdate;

  // std::vector< std::pair< HalfedgeGraph, bool > > meshes;
};

} // namespace FEVV

#include "Visualization/SimpleAdapterVisu.inl"
