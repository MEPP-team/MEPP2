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
class DialogCompressionValence1;
}
////////////////////////////////////////////////////////////////////////////////
namespace FEVV {

class DialogCompressionValence1 : public QDialog
{
  Q_OBJECT

public:
  explicit DialogCompressionValence1(QWidget *parent = 0);
  ~DialogCompressionValence1();

  void setCompressionValenceParams(const std::string &p3dFilePath,
                                   bool with_adaptative_quantization,
                                   int quantization_bits,
                                   int max_verticesdouble);
  void getCompressionValenceParams(std::string &p3dFilePath,
                                   bool &with_adaptative_quantization,
                                   int &quantization_bits,
                                   int &max_vertices);

private slots:
  void selectFilename();

private:
  Ui::DialogCompressionValence1 *ui;
};

} // namespace FEVV

////////////////////////////////////////////////////////////////////////////////
