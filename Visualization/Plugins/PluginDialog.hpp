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
