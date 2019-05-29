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

#if(FEVV_USE_QT5)
#include <QtWidgets/QApplication>
#else
#include <QtGui/QApplication>
#endif

#include <iostream>

namespace FEVV {

/**
 * \class SimpleApplication
 * \brief SimpleApplication is a specialization of QApplication.
 * This is useful if we want to catch exception.
 *
 * @see testViewer.cpp
 */
class SimpleApplication : public QApplication
{
public:
  /**
   * Constructor.
   */
  SimpleApplication(int &_argc, char **_argv) : QApplication(_argc, _argv)
  {
#ifdef DEBUG_VISU2
    std::cout << "*** this=" << this << "    entering " << __func__ << std::endl;
#endif

#ifdef DEBUG_VISU2
    std::cout << "*** this=" << this << "    leaving " << __func__ << std::endl;
#endif
  }

  ~SimpleApplication()
  {
#ifdef DEBUG_VISU2
    std::cout << "*** this=" << this << "    entering " << __func__ << std::endl;
#endif

#ifdef DEBUG_VISU2
    std::cout << "*** this=" << this << "    leaving " << __func__ << std::endl;
#endif
  }

private:
  /**
   * Default QApplication::notify() behavior.
   * but, we can catch exception.
   */
  bool notify(QObject *_receiver, QEvent *_event)
  {
#ifdef DEBUG_VISU2
    std::cout << "*** SimpleApplication notify event type=" << _event->type()
              << "  receiver=" << _receiver << std::endl;
#endif

    try
    {
      return QApplication::notify(_receiver, _event);
    }
    catch(std::exception &e)
    {
      // catch anything thrown within try block that derives from std::exception
      std::cerr << "std::exception was caught: " << e.what() << std::endl;

      return false;
    }
  }
};

} // namespace FEVV
