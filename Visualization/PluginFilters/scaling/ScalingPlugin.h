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
#include "Dialogs/DialogScaling1.h"

#ifndef Q_MOC_RUN // MT : very important to avoid the error : ' Parse error at
                  // "BOOST_JOIN" ' -> (qt4 pb with boost)
#include "Visualization/PluginFilters/BasePluginQt.h"
#include "Visualization/SimpleViewer.h"

#include "Visualization/SimpleWindow.h"

#include "FEVV/Filters/Generic/scaling.hpp" // A) include the header of the filter corresponding to your operation

#include "FEVV/Wrappings/properties.h"

#ifdef FEVV_USE_CGAL
#include "FEVV/Wrappings/properties_polyhedron_3.h"
#include "FEVV/Wrappings/properties_surface_mesh.h"
#include "FEVV/Wrappings/properties_linear_cell_complex.h"
#endif // FEVV_USE_CGAL
#ifdef FEVV_USE_OPENMESH
#include "FEVV/Wrappings/properties_openmesh.h"
#endif // FEVV_USE_OPENMESH
#ifdef FEVV_USE_AIF
#include "FEVV/Wrappings/properties_aif.h"
#endif // FEVV_USE_AIF
#endif

namespace FEVV {

class ScalingPlugin : public QObject,
                      public Generic_PluginInterface,
                      public BasePluginQt
{
  Q_OBJECT
  Q_INTERFACES(FEVV::Generic_PluginInterface)
#if(FEVV_USE_QT5) // see at the end of .cpp for QT4
  Q_PLUGIN_METADATA(IID "ScalingPlugin")
#endif

  /*public:
      using BasePlugin::apply;*/
public:
  ScalingPlugin() = default;
  ~ScalingPlugin() = default;

public:
  void init() override { init(true, 1.0, 1.0, 1.0); }

  void init(bool _forceCompute, double _x, double _y, double _z)
  {
    *value_forceCompute = _forceCompute;

    *value_x = _x;
    *value_y = _y;
    *value_z = _z;
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

    // window->setParam( "(Qt) Scale: X", value_x, "scaling_qt_p", this );
    // window->setParam( "(Qt) Scale: Y", value_y, "scaling_qt_p", this );
    // window->setParam( "(Qt) Scale: Z", value_z, "scaling_qt_p", this );
  }

  template< typename HalfedgeGraph >
  void scale(HalfedgeGraph *_mesh)
  {
    std::cout << "Asking to Scale mesh ! " << std::endl;

    auto pm = get(boost::vertex_point, *_mesh);

    Filters::calculate_scaling( // B) call the filter corresponding to your
                                // operation
        *_mesh,
        pm,
        *value_x,
        *value_y,
        *value_z);

    std::cout << "Scale mesh of " << *value_x << ";" << *value_y << ";"
              << *value_z << "." << std::endl;
  }

  template< typename HalfedgeGraph >
  void applyHG(BaseAdapterVisu *_adapter,
               HalfedgeGraph *_mesh,
               FEVV::PMapsContainer *pmaps_bag)
  {
    // get filter parameters from dialog window
    DialogScaling1 dial1;
    dial1.setScale(*value_x, *value_y, *value_z);
    if(dial1.exec() == QDialog::Accepted)
      dial1.getScale(*value_x, *value_y, *value_z);
    else
      return; // abort applying filter

    // apply filter
    scale(_mesh);

    // redraw mesh
    SimpleViewer *viewer =
        dynamic_cast< SimpleViewer * >(_adapter->getViewer());
    if(viewer)
      viewer->draw_or_redraw_mesh(_mesh, pmaps_bag, true, false);

    reset();

    viewer->frame();
  }

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

  void apply(BaseAdapterVisu *_adapter,
             CGALPointSet *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    applyHG< CGALPointSet >(_adapter, _mesh, pmaps_bag);
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

#ifdef FEVV_USE_PCL
  void apply(BaseAdapterVisu *_adapter,
             PCLPointCloud *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    applyHG< PCLPointCloud >(_adapter, _mesh, pmaps_bag);
  }
#endif


  QStringList Generic_plugins() const override
  {
    return QStringList() << "ScalingPlugin";
  }

  bool Generic_plugin(const QString &plugin) override
  {
    SimpleWindow *sw = static_cast< SimpleWindow * >(window);
      // dynamic_cast fails under OS X
    sw->onModificationParam("scaling_qt_p", this);
    sw->onApplyButton();

    return true;
  }

signals:
  void resetSignal();

protected:
  bool *value_forceCompute = new bool(false);

  double *value_x = new double(0.0);
  double *value_y = new double(0.0);
  double *value_z = new double(0.0);
};

} // namespace FEVV

