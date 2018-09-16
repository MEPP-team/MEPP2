#include "Visualization/Helpers/QtHelpers.h"

void
FEVV::Helpers::changeBackgroundColor(QWidget *_widget, const Color &_color)
{
  QPalette pal(_widget->palette());
  pal.setColor(
      QPalette::Background,
      QColor(_color.red(), _color.green(), _color.blue(), _color.alpha()));
  _widget->setAutoFillBackground(true);
  _widget->setPalette(pal);
}

void
FEVV::Helpers::changeTextColor(QWidget *_widget, const Color &_color)
{
  QPalette pal(_widget->palette());
  pal.setColor(
      QPalette::WindowText,
      QColor(_color.red(), _color.green(), _color.blue(), _color.alpha()));
  _widget->setAutoFillBackground(true);
  _widget->setPalette(pal);
}
