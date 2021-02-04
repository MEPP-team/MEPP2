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
#include "DialogCurvature1.h"
#include "ui_DialogCurvature1.h"

#include <QWhatsThis>
#include <QUrl>
#include <QDesktopServices>
////////////////////////////////////////////////////////////////////////////////
FEVV::DialogCurvature1::DialogCurvature1(QWidget *parent)
    : QDialog(parent), ui(new Ui::DialogCurvature1)
{
  ui->setupUi(this);

  // ---

  helpButton = new QPushButton("?");

  helpButton->setMaximumSize(32, 32);
  connect( helpButton, &QPushButton::clicked, []() {
    QWhatsThis::enterWhatsThisMode();
  } );

  ui->verticalLayout->addWidget(helpButton, 0, Qt::AlignRight);
}
////////////////////////////////////////////////////////////////////////////////
FEVV::DialogCurvature1::~DialogCurvature1()
{
  delete helpButton;

  delete ui;
}
////////////////////////////////////////////////////////////////////////////////
void
FEVV::DialogCurvature1::setCurvature(bool geod, double radius, bool Cmin_max, bool Dmin_max)
{
  Qt::CheckState state_geod =
      geod ? Qt::CheckState::Checked : Qt::CheckState::Unchecked;
  ui->checkBox_geod->setCheckState(state_geod);

  ui->lineEdit_radius->setText(QString::number(radius));

  Qt::CheckState state_Cmin_max =
      Cmin_max ? Qt::CheckState::Checked : Qt::CheckState::Unchecked;
  ui->checkBox_Cmin_max->setCheckState(state_Cmin_max);

  Qt::CheckState state_Dmin_max =
      Dmin_max ? Qt::CheckState::Checked : Qt::CheckState::Unchecked;
  ui->checkBox_Dmin_max->setCheckState(state_Dmin_max);
}
////////////////////////////////////////////////////////////////////////////////
void
FEVV::DialogCurvature1::getCurvature(bool &geod, double &radius, bool &Cmin_max, bool &Dmin_max)
{
  geod = ui->checkBox_geod->isChecked();

  radius = ui->lineEdit_radius->text().toDouble();

  Cmin_max = ui->checkBox_Cmin_max->isChecked();

  Dmin_max = ui->checkBox_Dmin_max->isChecked();
}
////////////////////////////////////////////////////////////////////////////////
bool FEVV::DialogCurvature1::event(QEvent *e)
{
  if (e->type() == QEvent::EnterWhatsThisMode)
  {
    QWhatsThis::leaveWhatsThisMode();

    QString link = "https://liris.cnrs.fr/mepp/doc/nightly/_filter_mesh_curvature.html";
    QDesktopServices::openUrl(QUrl(link));

    return true;
  }

  return QDialog::event(e);
}
////////////////////////////////////////////////////////////////////////////////
