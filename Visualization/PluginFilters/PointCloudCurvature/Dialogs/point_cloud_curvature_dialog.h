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
#include "Visualization/PluginFilters/BasePluginDialogQt.h"
////////////////////////////////////////////////////////////////////////////////
namespace Ui {
class PointCloudCurvatureDialog;
}
////////////////////////////////////////////////////////////////////////////////
namespace FEVV {

class PointCloudCurvatureDialog : public BasePluginDialogQt
{
  Q_OBJECT

public:
  explicit PointCloudCurvatureDialog(QWidget *parent = 0);
  ~PointCloudCurvatureDialog();

  void setParameters(unsigned int k, double radius, bool knn_search);
  void getParameters(unsigned int &k, double &radius, bool &knn_search);

private:
  Ui::PointCloudCurvatureDialog *ui;
};

} // namespace FEVV

////////////////////////////////////////////////////////////////////////////////
