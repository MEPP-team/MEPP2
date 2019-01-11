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

#include <QtPlugin>
#include <QMessageBox>

#ifndef Q_MOC_RUN // MT : very important to avoid the error : ' Parse error at
                  // "BOOST_JOIN" ' -> (qt4 pb with boost)
#include "Visualization/BaseWindowQt.h"
#endif

/*!
 * \brief Interfaces for plugins
 * These interfaces will be used for different plugins.
 */

namespace FEVV {

class Generic_PluginInterface
{
public:
  virtual ~Generic_PluginInterface() {}

  virtual QStringList Generic_plugins() const = 0;
  virtual bool Generic_plugin(const QString &plugin) = 0;

  virtual void init(BaseWindowQt *bwQt)
  {
    this->baseWindowQt = bwQt;
    // QMessageBox(QMessageBox::Information, "MEPP", "Generic_PluginInterface:
    // init.").exec();
  }

protected:
  BaseWindowQt *baseWindowQt; //!< the BaseWindowQt pointer
};

} // namespace FEVV

Q_DECLARE_INTERFACE(FEVV::Generic_PluginInterface,
                    "fr.liris.MEPP2.Generic_PluginInterface/1.0")
