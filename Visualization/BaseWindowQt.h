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

#include "Visualization/BaseWindow.h"
#include "Visualization/BaseAdapterVisuQt.h"

#include <QMainWindow>

namespace FEVV {

class BaseWindowQt : public QMainWindow, public BaseWindow
{
public:
  using AdapterQt = BaseAdapterVisuQt;

  using BaseWindow::isInit;
  using BaseWindow::isValid;

public:
  /**
   * Constructor.
   *
   * @note The main widget (with the scene) will be create and initialize if not
   * set (see other constructors).
   *
   * @param[in]   _parent   Pointer to the QWidget parent of this QWidget (used
   * by Qt) (Default value = 0).
   * @param[in]   _flags    Windows flags (used by Qt) (Default value = 0).
   */
  BaseWindowQt(QWidget *_parent = 0, Qt::WindowFlags _flags = 0)
      : BaseWindow(), QMainWindow(_parent, _flags)
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
  virtual ~BaseWindowQt()
  {
#ifdef DEBUG_VISU2
    std::cout << "*** this=" << this << "    entering " << __func__ << std::endl;
#endif

#ifdef DEBUG_VISU2
    std::cout << "*** this=" << this << "    leaving " << __func__ << std::endl;
#endif
  }
};

} // namespace FEVV
