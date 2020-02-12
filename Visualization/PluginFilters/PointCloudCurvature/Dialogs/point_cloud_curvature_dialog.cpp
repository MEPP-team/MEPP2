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
#include "point_cloud_curvature_dialog.h"
#include "ui_point_cloud_curvature_dialog.h"
////////////////////////////////////////////////////////////////////////////////
FEVV::PointCloudCurvatureDialog::PointCloudCurvatureDialog(QWidget *parent)
    : QDialog(parent), ui(new Ui::PointCloudCurvatureDialog)
{
  ui->setupUi(this);
}
////////////////////////////////////////////////////////////////////////////////
FEVV::PointCloudCurvatureDialog::~PointCloudCurvatureDialog() { delete ui; }
////////////////////////////////////////////////////////////////////////////////
void
FEVV::PointCloudCurvatureDialog::setParameters(unsigned int k, double radius, bool knn_search)
{
  ui->lineEdit_K->setText(QString::number(k));
  ui->lineEdit_Radius->setText(QString::number(radius));
  if(knn_search)
    ui->radioButton_knnsearch->setChecked(true);
  else
    ui->radioButton_radiussearch->setChecked(true);
}
////////////////////////////////////////////////////////////////////////////////
void
FEVV::PointCloudCurvatureDialog::getParameters(unsigned int &k, double &radius, bool &knn_search)
{
  k = ui->lineEdit_K->text().toUInt();
  radius = ui->lineEdit_Radius->text().toDouble();
  knn_search = ui->radioButton_knnsearch->isChecked();
}
////////////////////////////////////////////////////////////////////////////////
