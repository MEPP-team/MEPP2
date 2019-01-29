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

#include <QCheckBox>
#include <functional>

#ifndef Q_MOC_RUN // MT : very important to avoid the error : ' Parse error at
                  // "BOOST_JOIN" ' -> (qt4 pb with boost)
#include "Visualization/PluginFilters/BasePlugin.h"
#endif

namespace FEVV {
class SimpleCheckBox : public QCheckBox
{
  Q_OBJECT
public:
  std::string myPluginName;
  BasePlugin *myPlugin;

  bool *myBoolValue = nullptr;

  SimpleCheckBox(QWidget *parent = 0) : QCheckBox(parent) {}

  void setParams(bool *_value, std::string _pluginName, BasePlugin *_plugin)
  {
    this->myPluginName = _pluginName;
    this->myPlugin = _plugin;
    this->myBoolValue = _value;
  }
private slots:
  void modificationSlot()
  {
    emit modificationSignal(this->myPluginName, this->myPlugin);
    *myBoolValue = isChecked();
  }
  void resetSlot()
  {
    // std::cout << "resetSlot SimpleCheckBox: ";

    Qt::CheckState state =
        (*myBoolValue) ? Qt::CheckState::Checked : Qt::CheckState::Unchecked;
    setCheckState(state);

    // std::cout << *myBoolValue << std::endl;
  }
signals:
  void modificationSignal(std::string, BasePlugin *);
};

} // namespace FEVV
