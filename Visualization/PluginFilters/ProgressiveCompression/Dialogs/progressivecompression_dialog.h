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
#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <QDialog>
////////////////////////////////////////////////////////////////////////////////
namespace Ui {
class ProgressiveCompressionDialog;
}
////////////////////////////////////////////////////////////////////////////////
namespace FEVV {

class ProgressiveCompressionDialog : public QDialog
{
  Q_OBJECT

public:
  explicit ProgressiveCompressionDialog(QWidget *parent = 0);
  ~ProgressiveCompressionDialog();

  void getParameters(int &quantiz,
                     int &nb_batches,
                     int &min_vertices,
                     std::string &metric,
                     std::string &vkept,
                     std::string &predictor,
                     std::string &filepath,
                     std::string &batchstop);




private:
  Ui::ProgressiveCompressionDialog *ui;
};

} // namespace FEVV

////////////////////////////////////////////////////////////////////////////////
