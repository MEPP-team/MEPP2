#pragma once

#include <QLineEdit>
#include <functional>

#ifndef Q_MOC_RUN // MT : very important to avoid the error : ' Parse error at
                  // "BOOST_JOIN" ' -> (qt4 pb with boost)
#include "Visualization/PluginFilters/BasePlugin.h"
#endif

namespace FEVV {
class SimpleLineEdit : public QLineEdit
{
  Q_OBJECT
public:
  std::string myPluginName;
  BasePlugin *myPlugin;

  std::string *myStringValue = nullptr;
  int *myIntValue = nullptr;
  double *myDoubleValue = nullptr;
  float *myFloatValue = nullptr;

  SimpleLineEdit(QWidget *parent = 0) : QLineEdit(parent) {}

  void
  setParams(std::string *_value, std::string _pluginName, BasePlugin *_plugin)
  {
    this->myPluginName = _pluginName;
    this->myPlugin = _plugin;
    this->myStringValue = _value;
  }
  void setParams(int *_value, std::string _pluginName, BasePlugin *_plugin)
  {
    this->myPluginName = _pluginName;
    this->myPlugin = _plugin;
    this->myIntValue = _value;
  }
  void setParams(double *_value, std::string _pluginName, BasePlugin *_plugin)
  {
    this->myPluginName = _pluginName;
    this->myPlugin = _plugin;
    this->myDoubleValue = _value;
  }
  void setParams(float *_value, std::string _pluginName, BasePlugin *_plugin)
  {
    this->myPluginName = _pluginName;
    this->myPlugin = _plugin;
    this->myFloatValue = _value;
  }
private slots:
  void modificationSlot()
  {
    emit modificationSignal(this->myPluginName, this->myPlugin);
    if(myStringValue)
    {
      *myStringValue = text().toStdString();
    }
    else if(myIntValue)
    {
      *myIntValue = text().toInt();
    }
    else if(myDoubleValue)
    {
      *myDoubleValue = text().toDouble();
    }
    else if(myFloatValue)
    {
      *myFloatValue = text().toFloat();
    }
  }
  void resetSlot()
  {
    // std::cout << "resetSlot SimpleLineEdit: ";

    if(myStringValue)
    {
      setText(QString::fromUtf8((*myStringValue).c_str()));

      // std::cout << *myStringValue << std::endl;
    }
    else if(myIntValue)
    {
      setText(QString::fromUtf8(std::to_string(*myIntValue).c_str()));

      // std::cout << *myIntValue << std::endl;
    }
    else if(myDoubleValue)
    {
      setText(QString::fromUtf8(std::to_string(*myDoubleValue).c_str()));

      // std::cout << *myDoubleValue << std::endl;
    }
    else if(myFloatValue)
    {
      setText(QString::fromUtf8(std::to_string(*myFloatValue).c_str()));

      // std::cout << *myFloatValue << std::endl;
    }
  }
signals:
  void modificationSignal(std::string, BasePlugin *);
};

} // namespace FEVV
