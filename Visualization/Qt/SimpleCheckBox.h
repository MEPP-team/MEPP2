#pragma once

#if defined(SimpleCheckBox_RECURSES)
#error Recursive header files inclusion detected in SimpleCheckBox.h
#else // defined(SimpleCheckBox_RECURSES)
/** Prevents recursive inclusion of headers. */
#define SimpleCheckBox_RECURSES

#if !defined SimpleCheckBox_h
/** Prevents repeated inclusion of headers. */
#define SimpleCheckBox_h

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

#endif // !defined SimpleCheckBox_h

#undef SimpleCheckBox_RECURSES
#endif // else defined(SimpleCheckBox_RECURSES)
