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

#include <QDialog>

#include <QPushButton>
#include <QWhatsThis>
#include <QEvent>

#include <QUrl>
#include <QDesktopServices>

namespace FEVV {

/**
 * \brief  This class is intended to provide a help message
 *         to plugins dialogs.
 */
class BasePluginDialogQt : public QDialog
{
public:
  BasePluginDialogQt(QWidget *parent) : QDialog(parent)
  {
    helpButton = new QPushButton("?");
    helpButton->setMaximumSize(32, 32);

    // Qt5 only...
    /*connect( helpButton, &QPushButton::clicked, []() {
      QWhatsThis::enterWhatsThisMode();
    } );*/
  }

  ~BasePluginDialogQt() { delete helpButton; }

  void helpTriggered() { QDesktopServices::openUrl(QUrl(link)); }

  bool event(QEvent *e)
  {
    if (e->type() == QEvent::EnterWhatsThisMode)
    {
      QWhatsThis::leaveWhatsThisMode();

      helpTriggered();

      return true;
    }

    return QDialog::event(e);
  }

protected:
  QPushButton *helpButton;

  QString link;
};

} // namespace FEVV
