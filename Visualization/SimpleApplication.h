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
  SimpleApplication(int &_argc, char **_argv) : QApplication(_argc, _argv) {}

private:
  /**
   * Default QApplication::notify() behavior.
   * but, we can catch exception.
   */
  bool notify(QObject *_receiver, QEvent *_event)
  {
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
