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
