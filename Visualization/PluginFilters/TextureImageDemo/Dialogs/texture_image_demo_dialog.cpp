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
#include "texture_image_demo_dialog.h"
#include "ui_texture_image_demo_dialog.h"
////////////////////////////////////////////////////////////////////////////////
FEVV::TextureImageDemoDialog::TextureImageDemoDialog(QWidget *parent)
    : BasePluginDialogQt(parent), ui(new Ui::TextureImageDemoDialog)
{
  ui->setupUi(this);

  ui->verticalLayout->addWidget(helpButton, 0, Qt::AlignRight);
  QObject::connect( helpButton, SIGNAL(clicked(bool)), this, SLOT(onHelpTriggered()) );

  // ---

  link = "https://liris.cnrs.fr/mepp/doc/nightly"; // no help page for texture_image_demo_plugin for now...
}
////////////////////////////////////////////////////////////////////////////////
FEVV::TextureImageDemoDialog::~TextureImageDemoDialog() { delete ui; }
////////////////////////////////////////////////////////////////////////////////
void
FEVV::TextureImageDemoDialog::setParameters(/*double x, double y, double z*/)
{
  /*
  ui->lineEdit_X->setText(QString::number(x));
  ui->lineEdit_Y->setText(QString::number(y));
  ui->lineEdit_Z->setText(QString::number(z));
  */
}
////////////////////////////////////////////////////////////////////////////////
void
FEVV::TextureImageDemoDialog::getParameters(/*double &x, double &y, double &z*/)
{
  /*
  x = ui->lineEdit_X->text().toDouble();
  y = ui->lineEdit_Y->text().toDouble();
  z = ui->lineEdit_Z->text().toDouble();
  */
}
////////////////////////////////////////////////////////////////////////////////
