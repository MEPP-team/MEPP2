// Copyright (c) 2012-2020 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published 
// by the Free Software Foundation; either version 3 of the License, 
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#include "cmdm_dialog.h"
#include "ui_cmdm_dialog.h"
////////////////////////////////////////////////////////////////////////////////
FEVV::CMDMDialog::CMDMDialog(QWidget *parent)
    : BasePluginDialogQt(parent), ui(new Ui::CMDMDialog)
{
  ui->setupUi(this);

  ui->verticalLayout->addWidget(helpButton, 0, Qt::AlignRight);
  QObject::connect( helpButton, SIGNAL(clicked(bool)), this, SLOT(onHelpTriggered()) );

  // ---

  link = "https://liris.cnrs.fr/mepp/doc/nightly/_filter_c_m_d_m.html";
}
////////////////////////////////////////////////////////////////////////////////
FEVV::CMDMDialog::~CMDMDialog() { delete ui; }
////////////////////////////////////////////////////////////////////////////////
void
FEVV::CMDMDialog::getProcess(bool &one_two,
							 bool &two_one,
							 int &scales,
							 bool &colors)
{
	one_two = ui->radioButton_1to2->isChecked() || ui->radioButton_Sym->isChecked();
	two_one = ui->radioButton_2to1->isChecked() || ui->radioButton_Sym->isChecked();
	scales = ui->spinBox_scales->value();
    colors = ui->checkBox_color->isChecked();
}
////////////////////////////////////////////////////////////////////////////////
