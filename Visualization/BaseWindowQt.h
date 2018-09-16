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
  }

  /**
   * Destructor.
   */
  virtual ~BaseWindowQt() = default;
};

} // namespace FEVV
