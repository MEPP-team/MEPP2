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

#if(_MSC_VER >= 1400)
#ifndef _SCL_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif
#endif

#include "Visualization/Plugins/PluginInterface.h"

#include <QStringList>
#include "Dialogs/copy_graph_dialog.h"

#ifndef Q_MOC_RUN // MT : very important to avoid the error : ' Parse error at
                  // "BOOST_JOIN" ' -> (qt4 pb with boost)
#include "Visualization/PluginFilters/BasePluginQt.h"
#include "Visualization/SimpleViewer.h"

#include "Visualization/SimpleWindow.h"

// A) include the header of the filter corresponding to your operation
#include "FEVV/Wrappings/properties.h"
#include "FEVV/Filters/Generic/copy_graph.hpp"

#ifdef FEVV_USE_CGAL
#include "FEVV/Wrappings/properties_polyhedron_3.h"
#include "FEVV/Wrappings/properties_surface_mesh.h"
#include "FEVV/Wrappings/properties_linear_cell_complex.h"
#include "FEVV/Wrappings/properties_cgal_point_set.h"
#include "FEVV/Filters/CGAL/Point_set/copy_graph_cgal_point_set.hpp"
#endif // FEVV_USE_CGAL
#ifdef FEVV_USE_OPENMESH
#include "FEVV/Wrappings/properties_openmesh.h"
#endif // FEVV_USE_OPENMESH
#ifdef FEVV_USE_AIF
#include "FEVV/Wrappings/properties_aif.h"
#endif // FEVV_USE_AIF
#ifdef FEVV_USE_PCL
#include "FEVV/Wrappings/properties_pcl_point_cloud.h"
#include "FEVV/Filters/PCL/copy_graph_pcl.hpp"
#endif // FEVV_USE_PCL
#endif

namespace FEVV {

class CopyGraphPlugin : public QObject,
                         public Generic_PluginInterface,
                         public BasePluginQt
{
  Q_OBJECT
  Q_INTERFACES(FEVV::Generic_PluginInterface)
#if(FEVV_USE_QT5) // see at the end of .cpp for QT4
  Q_PLUGIN_METADATA(IID "CopyGraphPlugin")
#endif

  /*public:
      using BasePlugin::apply;*/
public:
  CopyGraphPlugin() = default;
  ~CopyGraphPlugin() = default;

public:
  void init() override
  {
    m_output_mesh_void = nullptr;
    m_output_pmaps_bag = nullptr;
  }

  void reset() override
  {
    init();

    emit resetSignal();
  }

  void addParameters(BaseWindow *_window) override
  {
    window = _window;
    if(window == nullptr || !window->isInit())
    {
      std::cerr << "BaseWindow is null or not initialized." << std::endl;
      return;
    }
  }

  template< typename HalfedgeGraphS,
            typename HalfedgeGraphT >
  void process(const HalfedgeGraphS       *_mesh,
               const FEVV::PMapsContainer *pmaps_bag)
  {
    std::cout << "Asking to CopyGraph mesh ! " << std::endl;

    // create output mesh and pmaps bag
    HalfedgeGraphT *output_mesh = new HalfedgeGraphT;
    m_output_pmaps_bag = new FEVV::PMapsContainer;

    // store output mesh for later display
    m_output_mesh_void = static_cast< void * >(output_mesh);

    // apply filter
    FEVV::Filters::copy_graph(
        *_mesh, *pmaps_bag, *output_mesh, *m_output_pmaps_bag);
  }


  template< typename HalfedgeGraphS,
            typename HalfedgeGraphT >
  void applyHG(BaseAdapterVisu      *_adapter,
               HalfedgeGraphS       *_mesh,
               FEVV::PMapsContainer *pmaps_bag)
  {
    // note: disable dialog box because no need for now
    //// get filter parameters from dialog window
    //CopyGraphDialog dialog;
    //dialog.setParameters(value_x, value_y, value_z);
    //if(dialog.exec() == QDialog::Accepted)
    //  dialog.getParameters(value_x, value_y, value_z);
    //else
    //  return; // abort applying filter

    // apply filter
    process< HalfedgeGraphS, HalfedgeGraphT >(_mesh, pmaps_bag);

    // redraw mesh
    SimpleViewer *viewer =
        dynamic_cast< SimpleViewer * >(_adapter->getViewer());
    if(viewer)
    {
      // draw output mesh
      auto output_mesh = static_cast< HalfedgeGraphT * >(m_output_mesh_void);
      if(output_mesh)
      {
        viewer->draw_or_redraw_mesh(output_mesh,
                                    m_output_pmaps_bag,
                                    false,
                                    false,
                                    "COPY");
      }
    }

    // comment next line to keep parameters values between calls
    reset();

    viewer->frame();
  }


#ifdef FEVV_USE_OPENMESH
  void apply(BaseAdapterVisu *_adapter,
             MeshOpenMesh *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    apply_choose_output< MeshOpenMesh >(_adapter, _mesh, pmaps_bag);
  }
#endif

#ifdef FEVV_USE_CGAL
  void apply(BaseAdapterVisu *_adapter,
             MeshLCC *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    apply_choose_output< MeshLCC >(_adapter, _mesh, pmaps_bag);
  }

  void apply(BaseAdapterVisu *_adapter,
             MeshSurface *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    apply_choose_output< MeshSurface >(_adapter, _mesh, pmaps_bag);
  }

  void apply(BaseAdapterVisu *_adapter,
             MeshPolyhedron *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    apply_choose_output< MeshPolyhedron >(_adapter, _mesh, pmaps_bag);
  }

  void apply(BaseAdapterVisu *_adapter,
             CGALPointSet *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    apply_choose_output< CGALPointSet >(_adapter, _mesh, pmaps_bag);
  }
#endif

#ifdef FEVV_USE_AIF
  void apply(BaseAdapterVisu *_adapter,
             MeshAIF *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    apply_choose_output< MeshAIF >(_adapter, _mesh, pmaps_bag);
  }
#endif

#ifdef FEVV_USE_PCL
  void apply(BaseAdapterVisu *_adapter,
             PCLPointCloud *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    apply_choose_output< PCLPointCloud >(_adapter, _mesh, pmaps_bag);
  }
#endif


  template< typename MeshT >
  void apply_choose_output(BaseAdapterVisu      *_adapter,
                           MeshT                *_mesh,
                           FEVV::PMapsContainer *pmaps_bag)
  {
    // ask the user for the datastructure
    std::string mesh_type = chooseDatastructureMsgBox();
    if(mesh_type == "NONE")
      return; // cancel pressed, aborting

    // apply plugin
    if(false)
    {
      // dummy case to be able to write all other cases as 'else if'
    }
#ifdef FEVV_USE_CGAL
    else if(mesh_type == "POLYHEDRON")
    {
      applyHG< MeshT, MeshPolyhedron >(_adapter, _mesh, pmaps_bag);
    }
    else if(mesh_type == "SURFACEMESH")
    {
      applyHG< MeshT, MeshSurface >(_adapter, _mesh, pmaps_bag);
    }
    else if(mesh_type == "LCC")
    {
      applyHG< MeshT, MeshLCC >(_adapter, _mesh, pmaps_bag);
    }
    else if(mesh_type == "CGALPOINTSET")
    {
      applyHG< MeshT, CGALPointSet >(_adapter, _mesh, pmaps_bag);
    }
#endif
#ifdef FEVV_USE_OPENMESH
    else if(mesh_type == "OPENMESH")
    {
      applyHG< MeshT, MeshOpenMesh >(_adapter, _mesh, pmaps_bag);
    }
#endif
#ifdef FEVV_USE_AIF
    else if(mesh_type == "AIF")
    {
      applyHG< MeshT, MeshAIF >(_adapter, _mesh, pmaps_bag);
    }
#endif
#ifdef FEVV_USE_PCL
    else if(mesh_type == "PCLPOINTCLOUD")
    {
      applyHG< MeshT, PCLPointCloud >(_adapter, _mesh, pmaps_bag);
    }
#endif
    else
    {
      QMessageBox::information(0,
                               "",
                               QObject::tr("Copy Graph is not yet compatible "
                                           "with the datastructure you have "
                                           "chosen!"));
    }
  }


  QStringList Generic_plugins() const override
  {
    return QStringList() << "CopyGraphPlugin";
  }


  bool Generic_plugin(const QString &/*plugin*/) override
  {
    SimpleWindow *sw = static_cast< SimpleWindow * >(window);
      // dynamic_cast fails under OS X
    sw->onModificationParam("copy_graph_qt_p", this);
      //TODO-elo-refactor-plugins
      // 1) the name of the function onModificationParam()
      //    is unrelated to its content!!!
      // 2) Generic_plugin() is called by SimpleWindow
      //    and calls back SimpleWindow functions here ;
      //    looks like a bad architecture
    sw->onApplyButton();

    return true;
  }

signals:
  void resetSignal();

protected:
  // filter output
  void *m_output_mesh_void;
  FEVV::PMapsContainer *m_output_pmaps_bag;
};

} // namespace FEVV

