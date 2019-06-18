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


#ifndef Q_MOC_RUN
#include "Visualization/PluginFilters/BasePlugin.h"
#endif

#include <QMessageBox>
#include <QPushButton>


namespace FEVV {

/**
 * \brief  This class is intended to provide some standard message boxes
 *         to all plugins.
 */
class BasePluginQt : public BasePlugin
{
public:

  // the apply(..., MeshT*, ...) functions below must be overridden
  // in derived class if the plugin supports the MeshT mesh type

#ifdef FEVV_USE_CGAL
  virtual void apply(BaseAdapterVisu *_adapter,
                     MeshPolyhedron *_mesh,
                     FEVV::PMapsContainer *pmaps_bag) override
  {
    QMessageBox::warning(
        0, "", QObject::tr("This filter is not compatible with Polyhedron_3!"));
  }


  virtual void apply(BaseAdapterVisu *_adapter,
                     MeshSurface *_mesh,
                     FEVV::PMapsContainer *pmaps_bag) override
  {
    QMessageBox::warning(
        0, "", QObject::tr("This filter is not compatible with Surface_mesh!"));
  }


  virtual void apply(BaseAdapterVisu *_adapter,
                     MeshLCC *_mesh,
                     FEVV::PMapsContainer *pmaps_bag) override
  {
    QMessageBox::warning(
        0, "", QObject::tr("This filter is not compatible with LCC!"));
  }


  virtual void apply(BaseAdapterVisu *_adapter,
                     CGALPointSet *_mesh,
                     FEVV::PMapsContainer *pmaps_bag) override
  {
    QMessageBox::warning(
        0, "", QObject::tr("This filter is not compatible with CGALPointSet!"));
  }
#endif


#ifdef FEVV_USE_OPENMESH
  virtual void apply(BaseAdapterVisu *_adapter,
                     MeshOpenMesh *_mesh,
                     FEVV::PMapsContainer *pmaps_bag) override
  {
    QMessageBox::warning(
        0, "", QObject::tr("This filter is not compatible with OpenMesh!"));
  }
#endif


#ifdef FEVV_USE_AIF
  virtual void apply(BaseAdapterVisu *_adapter,
                     MeshAIF *_mesh,
                     FEVV::PMapsContainer *pmaps_bag) override
  {
    QMessageBox::warning(
        0, "", QObject::tr("This filter is not compatible with AIF!"));
  }
#endif


#ifdef FEVV_USE_PCL
  virtual void apply(BaseAdapterVisu *_adapter,
                     PCLPointCloud *_mesh,
                     FEVV::PMapsContainer *pmaps_bag) override
  {
    QMessageBox::warning(
        0, "", QObject::tr("This filter is not compatible with PCLPointCloud!"));
  }
#endif


  // case where the plugin is applied when no mesh is opened
  virtual void apply(BaseAdapterVisu *_adapter,
                     void *_mesh_void,
                     FEVV::PMapsContainer *pmaps_bag) override
  {
    QMessageBox::warning(
        0, "", QObject::tr("To apply this filter, please first <b>open a mesh</b>!"));
  }


  // function duplicated here from SimpleWindow to get around link error
  //   Unresolved external symbol "public: static struct QMetaObject
  //     const FEVV::SimpleWindow::staticMetaObject"
  // on Visual C++ 2015
  // when SimpleWindow::chooseDatastructureMsgBox() is called from
  // within a plugin
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

    QPushButton *abortButton = msgbox.addButton(QMessageBox::Cancel);

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
};


} // namespace FEVV

