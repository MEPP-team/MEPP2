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


#include <string>

#include <QMessageBox>
#include <QPushButton>


namespace FEVV {

/**
 * Display the "choose your datastructure" message box.
 * Needed at least by SimpleWindow when opening a mesh and by
 * DecompressionValencePlugin when choosing the target datastructure.
 */
inline
std::string chooseDatastructureMsgBox(void)
{
  QMessageBox msgbox;
  msgbox.setWindowTitle("Datastructure");
  msgbox.setText("Choose a datastructure to store the mesh(es):");
  msgbox.setIcon(QMessageBox::Question);

  // here the Role is used to order the buttons
#ifdef FEVV_USE_CGAL
  QPushButton *polyhedron_button =
      msgbox.addButton("Polyhedron_3", QMessageBox::ResetRole);
  QPushButton *surfacemesh_button =
      msgbox.addButton("Surface_mesh", QMessageBox::ResetRole);
  QPushButton *lcc_button =
      msgbox.addButton("LCC", QMessageBox::ResetRole);
  QPushButton *cgalpointset_button =
      msgbox.addButton("CGALPointSet", QMessageBox::ResetRole);
#endif

#ifdef FEVV_USE_OPENMESH
  QPushButton *openmesh_button =
      msgbox.addButton("OpenMesh", QMessageBox::ResetRole);
#endif

#ifdef FEVV_USE_AIF
  QPushButton *aif_button = msgbox.addButton("AIF", QMessageBox::ResetRole);
#endif

#ifdef FEVV_USE_PCL
  QPushButton *pcl_button =
      msgbox.addButton("PCLPointCloud", QMessageBox::ResetRole);
#endif

  msgbox.addButton(QMessageBox::Cancel);

  msgbox.exec();
  std::string choice("NONE");

#ifdef FEVV_USE_CGAL
  if(msgbox.clickedButton() == polyhedron_button)
  {
    choice = "POLYHEDRON";
  }
  else if(msgbox.clickedButton() == surfacemesh_button)
  {
    choice = "SURFACEMESH";
  }
  else if(msgbox.clickedButton() == lcc_button)
  {
    choice = "LCC";
  }
  else if(msgbox.clickedButton() == cgalpointset_button)
  {
    choice = "CGALPOINTSET";
  }
#endif

#ifdef FEVV_USE_OPENMESH
  if(msgbox.clickedButton() == openmesh_button)
  {
    choice = "OPENMESH";
  }
#endif

#ifdef FEVV_USE_AIF
  if(msgbox.clickedButton() == aif_button)
  {
    choice = "AIF";
  }
#endif

#ifdef FEVV_USE_PCL
  if(msgbox.clickedButton() == pcl_button)
  {
    choice = "PCLPOINTCLOUD";
  }
#endif

  return choice;
}

} // namespace FEVV
