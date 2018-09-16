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
