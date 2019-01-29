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
#include "DialogMSDM21.h"
#include "ui_DialogMSDM21.h"
////////////////////////////////////////////////////////////////////////////////
FEVV::DialogMSDM21::DialogMSDM21(QWidget *parent)
    : QDialog(parent), ui(new Ui::DialogMSDM21)
{
  ui->setupUi(this);
}
////////////////////////////////////////////////////////////////////////////////
FEVV::DialogMSDM21::~DialogMSDM21() { delete ui; }
////////////////////////////////////////////////////////////////////////////////
void
FEVV::DialogMSDM21::getProcess(bool &one_two,
                               bool &two_one,
                               int &scales,
                               bool &colors)
{
  one_two =
      ui->radioButton_1to2->isChecked() || ui->radioButton_Sim->isChecked();
  two_one =
      ui->radioButton_2to1->isChecked() || ui->radioButton_Sim->isChecked();
  scales = ui->spinBox_scales->value();
  colors = ui->checkBox_color->isChecked();
}
////////////////////////////////////////////////////////////////////////////////
