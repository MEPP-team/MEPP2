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
#include "Base/Color.hpp"

#include <QThread>

namespace FEVV {

namespace Helpers {

/**
 * Change background color of the widget.
 *
 * @note It will also change the background color of the main widget.
 *
 * @param[in] _color the color to use.
 **/
void
changeBackgroundColor(QWidget *_widget, const Color &_color);

/**
 * Change text color of the widget.
 *
 * @note It will also change the text color of the main widget.
 *
 * @param[in] _color the color to use.
 **/
void
changeTextColor(QWidget *_widget, const Color &_color);

/*!
 * \class SleeperThread
 * \brief SleeperThread class.
 */
class SleeperThread : public QThread
{
public:
  /*!
   * \fn static void msleep(unsigned long msecs)
   * \brief Sleep function.
   *
   * \param msecs sleeping time in ms.
   */
  static void msleep(unsigned long msecs) { QThread::msleep(msecs); }
};

} // namespace Helpers

} // namespace FEVV
