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
class DialogDecompressionValence1;
}
////////////////////////////////////////////////////////////////////////////////
namespace FEVV {

class DialogDecompressionValence1 : public QDialog
{
  Q_OBJECT

public:
  explicit DialogDecompressionValence1(QWidget *parent = 0);
  ~DialogDecompressionValence1();

  void setDecompressionValenceParams(const std::string &p3dFilePath,
                                     int stop_level,
                                     bool write_info,
                                     bool write_intermediate_meshes,
                                     bool display_intermediate_meshes);
  void getDecompressionValenceParams(std::string &p3dFilePath,
                                     int &stop_level,
                                     bool &write_info,
                                     bool &write_intermediate_meshes,
                                     bool &display_intermediate_meshes);

private slots:
  void selectFilename();

private:
  Ui::DialogDecompressionValence1 *ui;
};

} // namespace FEVV

////////////////////////////////////////////////////////////////////////////////
