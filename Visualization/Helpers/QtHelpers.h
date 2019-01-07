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
