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

// #include "Visualization/BaseAdapterVisu.h"
// #include "Visualization/BaseWindow.h"

#include "Base/Color.hpp"
#include "Base/Assert.h"

// Generic mesh iterators
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <osg/PrimitiveSet>
#include "FEVV/Wrappings/Geometry_traits.h"

namespace FEVV {

enum RenderMethod {
  RENDER_POINTS = osg::PrimitiveSet::POINTS,
  RENDER_LINES = osg::PrimitiveSet::LINE_LOOP,
  RENDER_FILL = osg::PrimitiveSet::POLYGON
};

enum class RenderMode {
  RENDER_LEGACY = 0,
  RENDER_SHADERS_DIRECT_LIGHTING = 1,
  RENDER_SHADERS_INDIRECT_LIGHTING = 2
};

class BaseAdapterVisu;
class BaseWindow;

class BaseViewer
{
public:
  using Window = BaseWindow;
  using Adapter = BaseAdapterVisu;

public:
  /**
   * Constructor.
   */
  BaseViewer()
  {
#ifdef DEBUG_VISU2
    std::cout << "*** this=" << this << "    entering " << __func__ << std::endl;
#endif

    m_RenderMethod = RenderMethod::RENDER_FILL;
    m_RenderSuperimposedVertices = false;
    m_RenderSuperimposedVertices_Big = false;
    m_RenderSuperimposedEdges = false;

    m_UseVertexColor = false;
    m_UseFaceColor = false;
    m_UseTexture = false;

    m_Lighting = true;
    m_SmoothFlat_Shading = true;
    m_RenderMode = RenderMode::RENDER_LEGACY;

    m_ShowAxis = false;
    m_ShowGrid = false;
    m_Show_Vertex_Normals = false;
    m_ShowSelected = true;

    //

    m_ShowTranslateDragger = false;
    m_ShowRotateDragger = false;

    //

    m_redraw = false;
    m_recomputeNT_if_redraw = false;
    m_recreateOSGobj_if_redraw = true;
    m_step = 0.;

    //

    m_space_time = false;
    m_space_time_changeColorMode = true;

    m_time = false;

#ifdef DEBUG_VISU2
    std::cout << "*** this=" << this << "    leaving " << __func__ << std::endl;
#endif
  }

  virtual ~BaseViewer()
  {
#ifdef DEBUG_VISU2
    std::cout << "*** this=" << this << "    entering " << __func__ << std::endl;
#endif

#ifdef DEBUG_VISU2
    std::cout << "*** this=" << this << "    leaving " << __func__ << std::endl;
#endif
  }

  virtual void init() = 0;

  virtual bool isInit() const = 0;

  virtual bool isValid() const = 0;

  virtual bool isSelected() const { return bIsSelected; }

  virtual void setSelected(const bool _isSelected)
  {
    bIsSelected = _isSelected;
  }

  /**
   * attach a Adapter to the current Viewer.
   *
   * @note Must be called by adapter.init()
   */
  // virtual void attach( Adapter* _adapter ) = 0;
  void attach(Window *_window)
  {
    Assert::check(_window != nullptr,
                  "The given window is null.",
                  "BaseViewer::attach(Window)");
    myWindow = _window;
  }

  void attach(Adapter *_adapter)
  {
    Assert::check(_adapter != nullptr,
                  "The given adapter is null.",
                  "BaseViewer::attach(Adapter)");
    myAdapter = _adapter;
  }

  /**
   * Change background color of the scene.
   *
   * @param[in] _color the color to use.
   **/
  virtual bool changeBackgroundColor(const Color &_color) = 0;

  /**
   * Export the current view of the scene into a screenshot (PNG).
   *
   * @note Can not be const due to OSG.
   *
   * @param[in] _name name of the file to export (without extension).
   **/
  virtual bool saveScreenshot(const std::string &_name) = 0;

  Adapter *getAdapter() { return myAdapter; }
  Window *getWindow() { return myWindow; }

protected:
  Window *myWindow = nullptr;
  Adapter *myAdapter = nullptr;

  bool bIsInit = false;
  bool bIsSelected = false;

public:
  RenderMethod m_RenderMethod;
  bool m_RenderSuperimposedVertices;
  bool m_RenderSuperimposedVertices_Big;
  bool m_RenderSuperimposedEdges;

  bool m_UseVertexColor;
  bool m_UseFaceColor;
  bool m_UseTexture;

  bool m_Lighting;
  bool m_SmoothFlat_Shading;
  RenderMode m_RenderMode;

  bool m_ShowAxis;
  bool m_ShowGrid;
  bool m_Show_Vertex_Normals;
  bool m_ShowSelected;

  //

  bool m_ShowTranslateDragger;
  bool m_ShowRotateDragger;

  //

  bool m_redraw;
  bool m_recomputeNT_if_redraw;
  bool m_recreateOSGobj_if_redraw;
  float m_step;

  //

  bool m_space_time;
  bool m_space_time_changeColorMode;

  bool m_time;
};

} // namespace FEVV
