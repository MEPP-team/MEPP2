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
#include "Dialogs/point_cloud_curvature_dialog.h"

#ifndef Q_MOC_RUN // MT : very important to avoid the error : ' Parse error at
                  // "BOOST_JOIN" ' -> (qt4 pb with boost)
#include "Visualization/PluginFilters/BasePluginQt.h"
#include "Visualization/SimpleViewer.h"

#include "Visualization/SimpleWindow.h"

// A) include the header of the filter corresponding to your operation
#include "FEVV/Filters/Generic/PointCloud/point_cloud_curvature.hpp"

#include "FEVV/Wrappings/properties.h"

#ifdef FEVV_USE_CGAL
#include "FEVV/Wrappings/properties_polyhedron_3.h"
#include "FEVV/Wrappings/properties_surface_mesh.h"
#include "FEVV/Wrappings/properties_linear_cell_complex.h"
#include "FEVV/Wrappings/properties_cgal_point_set.h"
#endif // FEVV_USE_CGAL
#ifdef FEVV_USE_OPENMESH
#include "FEVV/Wrappings/properties_openmesh.h"
#endif // FEVV_USE_OPENMESH
#ifdef FEVV_USE_AIF
#include "FEVV/Wrappings/properties_aif.h"
#endif // FEVV_USE_AIF
#ifdef FEVV_USE_PCL
#include "FEVV/Wrappings/properties_pcl_point_cloud.h"
#endif // FEVV_USE_AIF
#endif // Q_MOC_RUN 

namespace FEVV {

class PointCloudCurvaturePlugin : public QObject,
                         public Generic_PluginInterface,
                         public BasePluginQt
{
  Q_OBJECT
  Q_INTERFACES(FEVV::Generic_PluginInterface)
#if(FEVV_USE_QT5) // see at the end of .cpp for QT4
  Q_PLUGIN_METADATA(IID "PointCloudCurvaturePlugin")
#endif

  /*public:
      using BasePlugin::apply;*/
public:
  PointCloudCurvaturePlugin() = default;
  ~PointCloudCurvaturePlugin() = default;

public:
  void init() override { m_k = 15; }

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

  template< typename HalfedgeGraph >
  void process(HalfedgeGraph *pc, FEVV::PMapsContainer *pmaps_bag)
  {
    std::cout << "Running filter PointCloudCurvature..." << std::endl;

    // retrieve geometry property map
    auto pm = get(boost::vertex_point, *pc);

    // create a vertex curvature property map
    auto v_curvm =
        FEVV::make_vertex_property_map< HalfedgeGraph, double >(*pc);

    // create a vertex color property map
    auto v_cm = FEVV::make_property_map(FEVV::vertex_color, *pc);
    FEVV::put_property_map(FEVV::vertex_color, *pc, *pmaps_bag, v_cm);

    // apply Point Cloud Curvature filter
    FEVV::Filters::point_cloud_curvature(*pc, pm, m_k, v_curvm, v_cm);

    std::cout << "Running filter PointCloudCurvature... done." << std::endl;
  }


  template< typename HalfedgeGraph >
  void applyHG(BaseAdapterVisu *_adapter,
               HalfedgeGraph *_mesh,
               FEVV::PMapsContainer *pmaps_bag)
  {
    // get filter parameters from dialog window
    PointCloudCurvatureDialog dialog;
    dialog.setParameters(m_k);
    if(dialog.exec() == QDialog::Accepted)
      dialog.getParameters(m_k);
    else
      return; // abort applying filter

    // apply filter
    process(_mesh, pmaps_bag);

    // redraw mesh
    SimpleViewer *viewer =
        dynamic_cast< SimpleViewer * >(_adapter->getViewer());
    if(viewer)
      viewer->draw_or_redraw_mesh(_mesh, pmaps_bag, true, false);

    // comment next line to keep parameters values between calls
    //reset();

    viewer->frame();
  }

#if 0 // Point cloud only
#ifdef FEVV_USE_OPENMESH
  void apply(BaseAdapterVisu *_adapter,
             MeshOpenMesh *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    applyHG< MeshOpenMesh >(_adapter, _mesh, pmaps_bag);
  }
#endif

#ifdef FEVV_USE_CGAL
  void apply(BaseAdapterVisu *_adapter,
             MeshLCC *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    applyHG< MeshLCC >(_adapter, _mesh, pmaps_bag);
  }

  void apply(BaseAdapterVisu *_adapter,
             MeshSurface *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    applyHG< MeshSurface >(_adapter, _mesh, pmaps_bag);
  }

  void apply(BaseAdapterVisu *_adapter,
             MeshPolyhedron *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    applyHG< MeshPolyhedron >(_adapter, _mesh, pmaps_bag);
  }
#endif

#ifdef FEVV_USE_AIF
  void apply(BaseAdapterVisu *_adapter,
             MeshAIF *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    applyHG< MeshAIF >(_adapter, _mesh, pmaps_bag);
  }
#endif
#endif // Point cloud only

#ifdef FEVV_USE_CGAL
  void apply(BaseAdapterVisu *_adapter,
             FEVV::CGALPointSet *pc,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    applyHG< FEVV::CGALPointSet >(_adapter, pc, pmaps_bag);
  }
#endif

#ifdef FEVV_USE_PCL
  void apply(BaseAdapterVisu *_adapter,
             FEVV::PCLPointCloud *pc,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    applyHG< FEVV::PCLPointCloud >(_adapter, pc, pmaps_bag);
  }
#endif

  QStringList Generic_plugins() const override
  {
    return QStringList() << "PointCloudCurvaturePlugin";
  }


  bool Generic_plugin(const QString &plugin) override
  {
    SimpleWindow *sw = static_cast< SimpleWindow * >(window);
      // dynamic_cast fails under OS X
    sw->onModificationParam("point_cloud_curvature_qt_p", this);
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
  // filter parameters
  unsigned int m_k;
};

} // namespace FEVV

