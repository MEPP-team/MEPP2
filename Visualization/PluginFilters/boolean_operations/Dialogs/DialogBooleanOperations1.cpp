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

#include "DialogBooleanOperations1.h"
#include "ui_DialogBooleanOperations1.h"
////////////////////////////////////////////////////////////////////////////////
FEVV::DialogBooleanOperations1::DialogBooleanOperations1(QWidget *parent)
    : QDialog(parent), ui(new Ui::DialogBooleanOperations1)
{
  ui->setupUi(this);
}
////////////////////////////////////////////////////////////////////////////////
FEVV::DialogBooleanOperations1::~DialogBooleanOperations1() { delete ui; }
////////////////////////////////////////////////////////////////////////////////
void
FEVV::DialogBooleanOperations1::setParameters(const std::string &operation)
{
  if( operation == "UNION")
    ui->radioButton_Union->setChecked(true);
  else if( operation == "INTER")
    ui->radioButton_Inter->setChecked(true);
  else
    ui->radioButton_Minus->setChecked(true);
}
////////////////////////////////////////////////////////////////////////////////
void
FEVV::DialogBooleanOperations1::getParameters(std::string &operation)
{
  if(ui->radioButton_Union->isChecked())
    operation = "UNION";
  else if(ui->radioButton_Inter->isChecked())
    operation = "INTER";
  else
    operation = "MINUS";
}
////////////////////////////////////////////////////////////////////////////////
