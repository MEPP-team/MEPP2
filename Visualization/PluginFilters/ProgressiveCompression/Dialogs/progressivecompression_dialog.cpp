// Copyright (c) 2012-2022 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#include "progressivecompression_dialog.h"
#include "ui_progressivecompression_dialog.h"
////////////////////////////////////////////////////////////////////////////////
FEVV::ProgressiveCompressionDialog::ProgressiveCompressionDialog(
    QWidget *parent)
    : QDialog(parent), ui(new Ui::ProgressiveCompressionDialog)
{
  ui->setupUi(this);
}
////////////////////////////////////////////////////////////////////////////////
FEVV::ProgressiveCompressionDialog::~ProgressiveCompressionDialog()
{
  delete ui;
}
////////////////////////////////////////////////////////////////////////////////
void
FEVV::ProgressiveCompressionDialog::getParameters(int &quantiz,
                                                    int &nb_max_batches,
                                                    int &min_vertices,
                                                    std::string &metric,
                                                    std::string &vkept,
                                                    std::string &predictor,
                                                    std::string &filepath,
                                                    std::string &batchstop)
{
  quantiz = ui->Quantiz->value();
  nb_max_batches = ui->Numberbatches->value();
  metric = ui->Metric->currentText().toStdString();
  vkept = ui->Operator->currentText().toStdString();
  predictor = ui->Predictor->currentText().toStdString();
  batchstop = ui->Batchstop->currentText().toStdString();
  filepath = ui->outputpath->text().toStdString();
  min_vertices = ui->Minimumvertices->value();

}
////////////////////////////////////////////////////////////////////////////////
