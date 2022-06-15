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
#include "Dialogs/weightedPCA_dialog.h"

#ifndef Q_MOC_RUN // MT : very important to avoid the error : ' Parse error at
                  // "BOOST_JOIN" ' -> (qt4 pb with boost)
#include "Visualization/PluginFilters/BasePluginQt.h"
#include "Visualization/SimpleViewer.h"

#include "Visualization/SimpleWindow.h"

// the header of the filter corresponding to the operation
#include "FEVV/Filters/Generic/PointCloud/WeightedPCANormals/point_cloud_normal_wpca.hpp"

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
#ifdef FEVV_USE_PCL
#include "FEVV/Wrappings/properties_pcl_point_cloud.h"
#endif // FEVV_USE_AIF
#endif // Q_MOC_RUN

namespace FEVV {

class WeightedPCAPlugin : public QObject,
                         public Generic_PluginInterface,
                         public BasePluginQt
{
  Q_OBJECT
  Q_INTERFACES(FEVV::Generic_PluginInterface)
#if(FEVV_USE_QT5) // see at the end of .cpp for QT4
  Q_PLUGIN_METADATA(IID "WeightedPCAPlugin")
#endif

  /*public:
      using BasePlugin::apply;*/
public:
  WeightedPCAPlugin() = default;
  ~WeightedPCAPlugin() = default;

public:
  void init() override { init(100, 0.01, 10000000); }

  void init(double _x, double _y, double _z)
  {
    n_neight = _x;
    noise = _y;
    curvature = _z;
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

  template< typename HalfedgeGraph >
  void process(HalfedgeGraph *cloud, FEVV::PMapsContainer *pmaps_bag, int n_neight, float noise, float curvature)
  {
    std::cout << "Running filter WeihgtedPCA normal computation..." << std::endl;

    std::cout << "create vertex-normal map" << std::endl;
    auto v_nm = make_property_map(FEVV::vertex_normal, *cloud);
    // store property map in property maps bag
    FEVV::put_property_map(FEVV::vertex_normal, *cloud, *pmaps_bag, v_nm);

    // retrieve point property map (aka geometry)
    auto pm = get(boost::vertex_point, *cloud);

    // apply filter
    noise = (float)noise/sqrt(3);
    noise = std::max(noise, 0.000001f);
    FEVV::Filters::compute_weighted_pca(*cloud, pm, v_nm, n_neight, noise, curvature);

    std::cout << "Running filter WeihgtedPCA normal computation... done." << std::endl;
  }

  template< typename HalfedgeGraph >
  void applyHG(BaseAdapterVisu *_adapter,
               HalfedgeGraph *_mesh,
               FEVV::PMapsContainer *pmaps_bag)
  {
    // get filter parameters from dialog window
    WeightedPCADialog dialog;
    dialog.setParameters(n_neight, noise, curvature);
    if(dialog.exec() == QDialog::Accepted)
      dialog.getParameters(n_neight, noise, curvature);
    else
      return; // abort applying filter

    // apply filter
    process(_mesh, pmaps_bag, n_neight, noise, curvature);

    SimpleViewer *viewer =
        dynamic_cast< SimpleViewer * >(_adapter->getViewer());
    if(viewer)
      viewer->draw_or_redraw_mesh(_mesh, pmaps_bag, true, false);

    // comment next line to keep parameters values between calls
    reset();

#if(FEVV_USE_QT5)
    // empty
#else
    viewer->frame(); // necessary or not ?
#endif
  }

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
    return QStringList() << "WeightedPCAPlugin";
  }

  bool Generic_plugin(const QString &/*plugin*/) override
  {
    SimpleWindow *sw = static_cast< SimpleWindow * >(window);
      // dynamic_cast fails under OS X
    sw->onModificationParam("weightedPCA_qt_p", this);
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
  double n_neight = 0;
  double noise = 0;
  double curvature = 0;
};

} // namespace FEVV

