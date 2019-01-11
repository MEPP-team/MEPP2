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

////////////////////////////////////////////////////////////////////////////////
#include <QDialog>
////////////////////////////////////////////////////////////////////////////////
namespace Ui {
class DialogProcessing1;
}
////////////////////////////////////////////////////////////////////////////////
namespace FEVV {

class DialogProcessing1 : public QDialog
{
  Q_OBJECT

public:
  explicit DialogProcessing1(QWidget *parent = 0);
  ~DialogProcessing1();

  void setProcess(double x, double y, double z);
  void getProcess(double &x, double &y, double &z);

private:
  Ui::DialogProcessing1 *ui;
};

} // namespace FEVV

////////////////////////////////////////////////////////////////////////////////
