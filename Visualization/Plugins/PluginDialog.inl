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
#include "Visualization/Plugins/PluginInterface.h"

#include <QPluginLoader>
#include <QStringList>
#include <QDir>

#include <QLabel>
#include <QGridLayout>
#include <QPushButton>
#include <QTreeWidget>
#include <QTreeWidgetItem>
#include <QHeaderView>

#include <iostream>

inline
FEVV::PluginDialog::PluginDialog(const QString &path,
                                 const QStringList &fileNames,
                                 QWidget *parent)
    : QDialog(parent), label(new QLabel), treeWidget(new QTreeWidget),
      okButton(new QPushButton(tr("OK")))
{
  treeWidget->setAlternatingRowColors(false);
  treeWidget->setSelectionMode(QAbstractItemView::NoSelection);
  treeWidget->setColumnCount(1);
  treeWidget->header()->hide();

  okButton->setDefault(true);

  connect(okButton, SIGNAL(clicked()), this, SLOT(close()));

  QGridLayout *mainLayout = new QGridLayout;
  mainLayout->setColumnStretch(0, 1);
  mainLayout->setColumnStretch(2, 1);
  mainLayout->addWidget(label, 0, 0, 1, 3);
  mainLayout->addWidget(treeWidget, 1, 0, 1, 3);
  mainLayout->addWidget(okButton, 2, 1);
  setLayout(mainLayout);

  interfaceIcon.addPixmap(style()->standardPixmap(QStyle::SP_DirOpenIcon),
                          QIcon::Normal,
                          QIcon::On);
  interfaceIcon.addPixmap(style()->standardPixmap(QStyle::SP_DirClosedIcon),
                          QIcon::Normal,
                          QIcon::Off);
  featureIcon.addPixmap(style()->standardPixmap(QStyle::SP_FileIcon));

  setWindowTitle(tr("Plugins information"));
  findPlugins(path, fileNames);
}

inline
void
FEVV::PluginDialog::findPlugins(const QString &path,
                                const QStringList &fileNames)
{
  label->setText(tr("Found the following plugins\n"
                    "(looked in %1):")
                     .arg(QDir::toNativeSeparators(path)));

  const QDir dir(path);
  foreach(QObject *plugin, QPluginLoader::staticInstances())
    populateTreeWidget(
        plugin,
        tr("%1 (Static Plugin)").arg(plugin->metaObject()->className()));

  foreach(QString fileName, fileNames)
  {
    std::cout << fileName.toStdString() << std::endl;
    QPluginLoader loader(dir.absoluteFilePath(fileName));
    QObject *plugin = loader.instance();
    if(plugin)
      populateTreeWidget(plugin, fileName);
  }
}

inline
void
FEVV::PluginDialog::populateTreeWidget(QObject *plugin, const QString &text)
{
  QTreeWidgetItem *pluginItem = new QTreeWidgetItem(treeWidget);
  pluginItem->setText(0, text);
  treeWidget->setItemExpanded(pluginItem, true);

  QFont boldFont = pluginItem->font(0);
  boldFont.setBold(true);
  pluginItem->setFont(0, boldFont);

  if(plugin)
  {
    Generic_PluginInterface *iGeneric_Filter =
        qobject_cast< Generic_PluginInterface * >(plugin);
    if(iGeneric_Filter)
      addItems(pluginItem,
               "Generic_PluginInterface",
               iGeneric_Filter->Generic_plugins());
  }
}

inline
void
FEVV::PluginDialog::addItems(QTreeWidgetItem *pluginItem,
                             const char *interfaceName,
                             const QStringList &features)
{
  QTreeWidgetItem *interfaceItem = new QTreeWidgetItem(pluginItem);
  interfaceItem->setText(0, interfaceName);
  interfaceItem->setIcon(0, interfaceIcon);

  foreach(QString feature, features)
  {
    if(feature.endsWith("..."))
      feature.chop(3);
    QTreeWidgetItem *featureItem = new QTreeWidgetItem(interfaceItem);
    featureItem->setText(0, feature);
    featureItem->setIcon(0, featureIcon);
  }
}
