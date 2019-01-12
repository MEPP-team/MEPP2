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
#include <QIcon>

QT_BEGIN_NAMESPACE
class QLabel;
class QPushButton;
class QStringList;
class QTreeWidget;
class QTreeWidgetItem;
QT_END_NAMESPACE

namespace FEVV {

class PluginDialog : public QDialog
{
  Q_OBJECT

public:
  PluginDialog(const QString &path,
               const QStringList &fileNames,
               QWidget *parent = 0);

private:
  void findPlugins(const QString &path, const QStringList &fileNames);
  void populateTreeWidget(QObject *plugin, const QString &text);
  void addItems(QTreeWidgetItem *pluginItem,
                const char *interfaceName,
                const QStringList &features);

  QLabel *label;
  QTreeWidget *treeWidget;
  QPushButton *okButton;
  QIcon interfaceIcon;
  QIcon featureIcon;
};

} // namespace FEVV


#ifndef Q_MOC_RUN
// implementation
#include "PluginDialog.inl"
#endif // Q_MOC_RUN
