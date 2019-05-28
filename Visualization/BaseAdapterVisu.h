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

#include "Base/Assert.h"

// #include "Visualization/BaseWindow.h"
#include "Visualization/PluginFilters/BasePlugin.h"
#include "Visualization/BaseViewer.h"

#include <functional>

namespace FEVV {

class BaseViewer;
class BasePlugin;

class BaseAdapterVisu
{
public:
  using Viewer = BaseViewer;
  using Plugin = BasePlugin;

public:
  /**
   * Constructor.
   *
   * @note The OpenSceneGraph viewer will be created and initialized if not set
   * (see other constructors).
   */
  BaseAdapterVisu()
  {
#ifdef DEBUG_VISU2
    std::cout << "*** this=" << this << "    entering " << __func__ << std::endl;
#endif

#ifdef DEBUG_VISU2
    std::cout << "*** this=" << this << "    leaving " << __func__ << std::endl;
#endif
  }

  virtual ~BaseAdapterVisu()
  {
#ifdef DEBUG_VISU2
    std::cout << "*** this=" << this << "    entering " << __func__ << std::endl;
#endif

#ifdef DEBUG_VISU2
    std::cout << "*** this=" << this << "    leaving " << __func__ << std::endl;
#endif
  }

  /**
   * attach a Viewer to the current Widget.
   */
  virtual void attach(Viewer *_viewer)
  {
    Assert::check(_viewer != nullptr,
                  "The given viewer is null.",
                  "BaseAdapterVisu::attach(Viewer)");
    myViewer = _viewer;
    _viewer->attach(this);
  }

  /**
   * attach a Window to the current Widget.
   *
   * @note Must be called by window.init()
   */
  // virtual void attach( Window* _window ) = 0;

  /**
   * initialize the widget.
   *
   * @note Must attach himself to the viewer.
   */
  virtual void init() = 0;

  virtual bool isInit() const = 0;

  virtual bool isValid() const = 0;

  virtual bool isSelected() const
  {
    if(myViewer != nullptr)
    {
      return myViewer->isSelected();
    }
    else
    {
      Assert::check(myViewer == nullptr,
                    "myViewer is not attached.",
                    "SimpleWidget::setSelected");
      return false;
    }
  }

  virtual void setSelected(const bool _isSelected)
  {
    if(myViewer != nullptr)
    {
      return myViewer->setSelected(_isSelected);
    }
    else
    {
      Assert::check(myViewer == nullptr,
                    "myViewer is not attached.",
                    "SimpleWidget::setSelected");
    }
  }

  // virtual bool hasWindow() const
  // {
  //     return (myWindow != nullptr);
  // }

  /**
   * Get the viewer.
   *
   * @return the viewer.
   */
  virtual Viewer *getViewer()
  {
    Assert::check(myViewer != nullptr,
                  "No viewer attached.",
                  "BaseAdapterVisu::getViewer()");
    return myViewer;
  }

  // virtual Viewer* getViewer( int _position ) = 0;

  // std::vector< Viewer* > getViewer()
  // {
  //     return myViewers;
  // }

  // virtual std::vector< Viewer* > getSelectedViewers() = 0;
  // {
  //     std::vector< Viewer* > v_result;
  //     for( unsigned int iViewer = 0; iViewer < myViewers.size(); ++iViewer )
  //     {
  //         if( myViewers[ iViewer ]->isSelected() )
  //         {
  //             v_result.push_back( myViewers[ iViewer ] );
  //         }
  //     }
  //     return v_result;
  // }

  // unsigned int getNbViewer()
  // {
  //   return myViewers.size();
  // }

  /**
   * Get the window.
   *
   * @return the window.
   */
  // virtual Window* getWindow() = 0;

  virtual void apply(Plugin *myPlugin) = 0;

protected:
  Viewer *myViewer = nullptr;
  // Window* myWindow = nullptr;
  bool bIsInit = false;
  // bool bIsSelected = false;
};

} // namespace FEVV
