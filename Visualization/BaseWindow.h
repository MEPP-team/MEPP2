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

#include "Visualization/PluginFilters/BasePlugin.h"
#include "Visualization/BaseAdapterVisu.h"
// #include "Visualization/BaseViewer.h"

namespace FEVV {

class BaseAdapterVisu;
class BasePlugin;

class BaseWindow
{
public:
  using Adapter = BaseAdapterVisu;
  using Plugin = BasePlugin;

public:
  /**
   * Constructor.
   *
   * @param[in]   _parent   Pointer to the QWidget parent of this QWidget (used
   * by Qt) (Default value = 0).
   * @param[in]   _flags    Windows flags (used by Qt) (Default value = 0).
   */
  BaseWindow()
  {
#ifdef DEBUG_VISU2
    std::cout << "*** this=" << this << "    entering " << __func__ << std::endl;
#endif

#ifdef DEBUG_VISU2
    std::cout << "*** this=" << this << "    leaving " << __func__ << std::endl;
#endif
  }

  /**
   * Destructor.
   */
  virtual ~BaseWindow()
  {
#ifdef DEBUG_VISU2
    std::cout << "*** this=" << this << "    entering " << __func__ << std::endl;
#endif

#ifdef DEBUG_VISU2
    std::cout << "*** this=" << this << "    leaving " << __func__ << std::endl;
#endif
  }

  /**
   * Attach a Visualization Adapter to the current Window.
   *
   * @param[in]   _adapter   Pointer to a visualization adapter (with the scene,
   * see BaseAdapterVisu).
   */
  virtual void attach(Adapter *_adapter) = 0;

  virtual void attachPlugin(Plugin *_plugin) = 0;

  virtual void setParam(std::string _name,
                        int *_value,
                        std::string _pluginName,
                        Plugin *_plugin) = 0;

  virtual void setParam(std::string _name,
                        double *_value,
                        std::string _pluginName,
                        Plugin *_plugin) = 0;

  virtual void setParam(std::string _name,
                        float *_value,
                        std::string _pluginName,
                        Plugin *_plugin) = 0;

  virtual void setParam(std::string _name,
                        bool *_value,
                        std::string _pluginName,
                        Plugin *_plugin) = 0;

  virtual void setParam(std::string _name,
                        std::string *_value,
                        std::string _pluginName,
                        Plugin *_plugin) = 0;

  /**
   * initialize the window.
   */
  virtual void init() = 0;

  virtual bool isInit() const { return bIsInit; }

  virtual bool isValid() const { return bIsInit && (adapters.size() != 0); }

  virtual void notify() = 0;

  /**
   * Get all visualization adapters.
   * @return the pointer of the vector of adapters.
   */
  virtual std::vector< Adapter * > *getAdapters() { return &adapters; }

  virtual Adapter *getAdapter(unsigned int _pos)
  {
    Assert::check(_pos < adapters.size(),
                  "No adapter attached at this position.",
                  "BaseWindow::getAdapter(int)");
    return adapters[_pos];
  }

  virtual std::vector< Adapter * > getSelectedAdapters() = 0;
  virtual std::vector< Adapter::Viewer * > getSelectedViewers() = 0;

protected:
  std::vector< Adapter * > adapters;
  std::map< std::string, Plugin * > stackPlugins;
  bool bIsInit = false;
};

} // namespace FEVV
