#pragma once

#include "Visualization/BaseAdapterVisu.h"
#include "Visualization/BaseViewerOSG.h"

#include <QWidget>
#include <QGridLayout>

namespace FEVV {

class BaseAdapterVisuQt : public QWidget, public BaseAdapterVisu
{
public:
  using ViewerOSG = BaseViewerOSG;

  using BaseAdapterVisu::isSelected;
  using BaseAdapterVisu::setSelected;

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
  BaseAdapterVisuQt(QWidget *_parent = 0, Qt::WindowFlags _f = 0)
      : BaseAdapterVisu(), QWidget(_parent, _f)
  {
    layout = new QGridLayout(this);
    layout->setObjectName(QString::fromUtf8("WidgetLayout"));
  }

  virtual ~BaseAdapterVisuQt() { delete layout; }

  void attach(Viewer *_viewer) override
  {
    ViewerOSG *viewer = dynamic_cast< ViewerOSG * >(_viewer);

    if(viewer)
    {
      attach(viewer);
    }
    else
    {

      Assert::check(false,
                    "is not implemented. See attach(ViewerOSG) instead",
                    "BaseAdapterVisuQt::attach(Viewer)");
    }
  }

  virtual void attach(ViewerOSG *_viewer)
  {
    Assert::check(_viewer != nullptr,
                  "The given viewer is null",
                  "BaseAdapterVisuQt::attach(ViewerOSG)");
    myViewer = _viewer;
    _viewer->attach(this);
  }

protected:
  QGridLayout *layout;
};

} // namespace FEVV
