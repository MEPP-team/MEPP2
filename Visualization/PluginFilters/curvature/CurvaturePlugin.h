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
#include "Dialogs/DialogCurvature1.h"

#ifndef Q_MOC_RUN // MT : very important to avoid the error : ' Parse error at
                  // "BOOST_JOIN" ' -> (qt4 pb with boost)
#include "Visualization/PluginFilters/BasePluginQt.h"
#include "Visualization/SimpleViewer.h"

#include "Visualization/SimpleWindow.h"

#include "FEVV/Filters/Generic/Manifold/Curvature/curvature.hpp" // A) include the header of the filter corresponding to your operation
#include "FEVV/Filters/Generic/calculate_face_normals.hpp"
#include "FEVV/Filters/Generic/AABB/get_max_bb_size.hpp"

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
#endif // Q_MOC_RUN

namespace FEVV {

class CurvaturePlugin : public QObject,
                        public Generic_PluginInterface,
                        public BasePluginQt
{
  Q_OBJECT
  Q_INTERFACES(FEVV::Generic_PluginInterface)
#if(FEVV_USE_QT5) // see at the end of .cpp for QT4
  Q_PLUGIN_METADATA(IID "CurvaturePlugin")
#endif

  /*public:
      using BasePlugin::apply;*/
public:
  CurvaturePlugin() { initLUT(); }

  ~CurvaturePlugin() = default;

public:
  void init() override { init(true, false, 0.001, 2, false, false); }

  void init(bool _forceCompute,
            bool _isGeod,
            double _radius,
            int _colorField,
            bool _displayMinDirections,
            bool _displayMaxDirections)
  {
    *value_forceCompute = _forceCompute;

    *value_isGeod = _isGeod;
    *value_radius = _radius;

    *value_colorField = _colorField; // 0 : rien, 1 : min, 2 : max

    *value_displayMinDirections = _displayMinDirections;
    *value_displayMaxDirections = _displayMaxDirections;
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

    // window->setParam( "(Qt) Curvature: geodesic", value_isGeod,
    // "curvature_qt_p", this ); window->setParam( "(Qt) Curvature: radius",
    // value_radius, "curvature_qt_p", this );
  }

  void initLUT()
  {
    int i = 0;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.515600;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.531300;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.546900;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.562500;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.578100;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.593800;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.609400;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.625000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.640600;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.656300;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.671900;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.687500;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.703100;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.718800;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.734400;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.750000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.765600;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.781300;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.796900;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.812500;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.828100;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.843800;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.859400;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.875000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.890600;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.906300;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.921900;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.937500;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.953100;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.968800;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.984400;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.015600;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.031300;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.046900;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.062500;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.078100;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.093800;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.109400;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.125000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.140600;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.156300;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.171900;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.187500;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.203100;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.218800;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.234400;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.250000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.265600;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.281300;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.296900;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.312500;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.328100;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.343800;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.359400;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.375000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.390600;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.406300;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.421900;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.437500;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.453100;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.468800;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.484400;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.500000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.515600;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.531300;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.546900;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.562500;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.578100;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.593800;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.609400;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.625000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.640600;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.656300;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.671900;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.687500;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.703100;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.718800;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.734400;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.750000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.765600;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.781300;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.796900;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.812500;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.828100;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.843800;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.859400;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.875000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.890600;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.906300;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.921900;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.937500;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.953100;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.968800;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.984400;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.015600;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.031300;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.984400;
    LUT_CourbureClust[i++] = 0.046900;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.968800;
    LUT_CourbureClust[i++] = 0.062500;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.953100;
    LUT_CourbureClust[i++] = 0.078100;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.937500;
    LUT_CourbureClust[i++] = 0.093800;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.921900;
    LUT_CourbureClust[i++] = 0.109400;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.906300;
    LUT_CourbureClust[i++] = 0.125000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.890600;
    LUT_CourbureClust[i++] = 0.140600;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.875000;
    LUT_CourbureClust[i++] = 0.156300;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.859400;
    LUT_CourbureClust[i++] = 0.171900;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.843800;
    LUT_CourbureClust[i++] = 0.187500;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.828100;
    LUT_CourbureClust[i++] = 0.203100;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.812500;
    LUT_CourbureClust[i++] = 0.218800;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.796900;
    LUT_CourbureClust[i++] = 0.234400;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.781300;
    LUT_CourbureClust[i++] = 0.250000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.765600;
    LUT_CourbureClust[i++] = 0.265600;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.750000;
    LUT_CourbureClust[i++] = 0.281300;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.734400;
    LUT_CourbureClust[i++] = 0.296900;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.718800;
    LUT_CourbureClust[i++] = 0.312500;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.703100;
    LUT_CourbureClust[i++] = 0.328100;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.687500;
    LUT_CourbureClust[i++] = 0.343800;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.671900;
    LUT_CourbureClust[i++] = 0.359400;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.656300;
    LUT_CourbureClust[i++] = 0.375000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.640600;
    LUT_CourbureClust[i++] = 0.390600;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.625000;
    LUT_CourbureClust[i++] = 0.406300;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.609400;
    LUT_CourbureClust[i++] = 0.421900;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.593800;
    LUT_CourbureClust[i++] = 0.437500;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.578100;
    LUT_CourbureClust[i++] = 0.453100;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.562500;
    LUT_CourbureClust[i++] = 0.468800;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.546900;
    LUT_CourbureClust[i++] = 0.484400;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.531300;
    LUT_CourbureClust[i++] = 0.500000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.515600;
    LUT_CourbureClust[i++] = 0.515600;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.500000;
    LUT_CourbureClust[i++] = 0.531300;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.484400;
    LUT_CourbureClust[i++] = 0.546900;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.468800;
    LUT_CourbureClust[i++] = 0.562500;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.453100;
    LUT_CourbureClust[i++] = 0.578100;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.437500;
    LUT_CourbureClust[i++] = 0.593800;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.421900;
    LUT_CourbureClust[i++] = 0.609400;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.406300;
    LUT_CourbureClust[i++] = 0.625000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.390600;
    LUT_CourbureClust[i++] = 0.640600;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.375000;
    LUT_CourbureClust[i++] = 0.656300;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.359400;
    LUT_CourbureClust[i++] = 0.671900;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.343800;
    LUT_CourbureClust[i++] = 0.687500;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.328100;
    LUT_CourbureClust[i++] = 0.703100;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.312500;
    LUT_CourbureClust[i++] = 0.718800;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.296900;
    LUT_CourbureClust[i++] = 0.734400;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.281300;
    LUT_CourbureClust[i++] = 0.750000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.265600;
    LUT_CourbureClust[i++] = 0.765600;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.250000;
    LUT_CourbureClust[i++] = 0.781300;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.234400;
    LUT_CourbureClust[i++] = 0.796900;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.218800;
    LUT_CourbureClust[i++] = 0.812500;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.203100;
    LUT_CourbureClust[i++] = 0.828100;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.187500;
    LUT_CourbureClust[i++] = 0.843800;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.171900;
    LUT_CourbureClust[i++] = 0.859400;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.156300;
    LUT_CourbureClust[i++] = 0.875000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.140600;
    LUT_CourbureClust[i++] = 0.890600;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.125000;
    LUT_CourbureClust[i++] = 0.906300;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.109400;
    LUT_CourbureClust[i++] = 0.921900;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.093800;
    LUT_CourbureClust[i++] = 0.937500;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.078100;
    LUT_CourbureClust[i++] = 0.953100;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.062500;
    LUT_CourbureClust[i++] = 0.968800;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.046900;
    LUT_CourbureClust[i++] = 0.984400;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.031300;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.015600;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.984400;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.968800;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.953100;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.937500;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.921900;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.906300;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.890600;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.875000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.859400;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.843800;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.828100;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.812500;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.796900;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.781300;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.765600;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.750000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.734400;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.718800;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.703100;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.687500;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.671900;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.656300;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.640600;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.625000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.609400;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.593800;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.578100;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.562500;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.546900;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.531300;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.515600;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.500000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.484400;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.468800;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.453100;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.437500;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.421900;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.406300;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.390600;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.375000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.359400;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.343800;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.328100;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.312500;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.296900;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.281300;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.265600;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.250000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.234400;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.218800;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.203100;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.187500;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.171900;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.156300;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.140600;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.125000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.109400;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.093800;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.078100;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.062500;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.046900;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.031300;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.015600;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 1.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.984400;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.968800;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.953100;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.937500;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.921900;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.906300;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.890600;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.875000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.859400;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.843800;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.828100;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.812500;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.796900;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.781300;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.765600;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.750000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.734400;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.718800;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.703100;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.687500;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.671900;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.656300;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.640600;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.625000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.609400;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.593800;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.578100;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.562500;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.546900;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.531300;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.515600;
    LUT_CourbureClust[i++] = 0.000000;
    LUT_CourbureClust[i++] = 0.000000;
  }

  template< typename HalfedgeGraph,
            typename VertexCurvatureMap,
            typename VertexColorMap,
            typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
  void constructColorMap(const HalfedgeGraph &g,
                         VertexCurvatureMap &v_cm,
                         VertexColorMap &v_colorm,
                         double MinNrmMinCurvature,
                         double MaxNrmMinCurvature,
                         double MinNrmMaxCurvature,
                         double MaxNrmMaxCurvature,
                         int ColorField) // 0 : rien, 1 : min, 2 : max
  {
    typedef boost::graph_traits< HalfedgeGraph > GraphTraits;
    typedef typename GraphTraits::vertex_iterator vertex_iterator;

    using Vector = typename GeometryTraits::Vector;

    double R;
    int indiceLut;

    auto iterator_pair =
        vertices(g); // vertices() returns a vertex_iterator pair
    vertex_iterator vi = iterator_pair.first;
    vertex_iterator vi_end = iterator_pair.second;
    for(; vi != vi_end; ++vi)
    {
      if(ColorField == 1)
        R = (v_cm[*vi].KminCurv - MinNrmMinCurvature) /
            (MaxNrmMinCurvature - MinNrmMinCurvature) * 255;
      else if(ColorField == 2)
        R = (v_cm[*vi].KmaxCurv - MinNrmMaxCurvature) /
            (MaxNrmMaxCurvature - MinNrmMaxCurvature) * 255;
      else
        R = 1;

      if(R > 255)
        R = 255;
      indiceLut = floor(R);

      v_colorm[*vi] = Vector(LUT_CourbureClust[3 * indiceLut],
                             LUT_CourbureClust[3 * indiceLut + 1],
                             LUT_CourbureClust[3 * indiceLut + 2]);
    }
  }

  template< typename HalfedgeGraph, typename VertexCurvatureMap, typename GuiPropertiesMap >
  void curvature(HalfedgeGraph *_mesh,
                 GuiPropertiesMap &m_gpm,
                 VertexCurvatureMap &v_cm,
                 double &MinNrmMinCurvature,
                 double &MaxNrmMinCurvature,
                 double &MinNrmMaxCurvature,
                 double &MaxNrmMaxCurvature)
  {
    std::cout << "Asking to Curvature mesh ! " << std::endl;

    auto pm = get(boost::vertex_point, *_mesh);

    double max_bb_size = Filters::get_max_bb_size(*_mesh, pm);

    auto f_nm = make_property_map(FEVV::face_normal, *_mesh);

    Filters::calculate_face_normals(*_mesh, pm, f_nm);

    Filters::calculate_curvature( // B) call the filter corresponding to your
                                  // operation
        *_mesh,
        v_cm,
        pm,
        f_nm,
        *value_isGeod,
        *value_radius * max_bb_size,
        MinNrmMinCurvature,
        MaxNrmMinCurvature,
        MinNrmMaxCurvature,
        MaxNrmMaxCurvature); // minimum and maximum values of the minimum and
                             // maximum curvature fields (usefull for color
                             // rendering)

    // change display mode to show filter result
    auto gui_props = get(m_gpm, 0);
    //gui_props.display_mode = FEVV::Types::DisplayMode::VERTEX_COLOR; // not used yet...
    //gui_props.is_visible = true; // not necessary because true by default...
    put(m_gpm, 0, gui_props);

    std::cout << "Curvature mesh, isGeod: " << *value_isGeod
              << " - radius: " << *value_radius << "." << std::endl;
  }

  //#define DBG_applyHG

  template< typename HalfedgeGraph >
  void applyHG(BaseAdapterVisu *_adapter,
               HalfedgeGraph *_mesh,
               FEVV::PMapsContainer *pmaps_bag)
  {
    // get filter parameters from dialog window
    DialogCurvature1 dial1;
    dial1.setCurvature(*value_isGeod, *value_radius);
    if(dial1.exec() == QDialog::Accepted)
      dial1.getCurvature(*value_isGeod, *value_radius);
    else
      return; // abort applying filter

    using VertexCurvatureMapHG =
        FEVV::Vertex_pmap< HalfedgeGraph, Filters::v_Curv< HalfedgeGraph > >;

    static std::map<
        std::pair< BaseAdapterVisu *, HalfedgeGraph * >,
        std::tuple< VertexCurvatureMapHG, std::array< double, 4 > > >
        map_v_cmHG; // static to allow to not recompute

    // create property map to store vertex color
    auto v_colormHG = make_property_map(FEVV::vertex_color, *_mesh);
    put_property_map(FEVV::vertex_color, *_mesh, *pmaps_bag, v_colormHG);

    // sample/nothing to do with curvature - create property map to store edge
    // color
    using VectorHG = typename Geometry_traits< HalfedgeGraph >::Vector;

    auto e_colormHG = make_property_map(FEVV::edge_color, *_mesh);
    put_property_map(FEVV::edge_color, *_mesh, *pmaps_bag, e_colormHG);

    typedef boost::graph_traits< HalfedgeGraph > GraphTraits;
    typedef typename GraphTraits::edge_iterator edge_iterator;
    edge_iterator eb, ee;
    for(boost::tie(eb, ee) = edges(*_mesh); eb != ee; ++eb)
      put(e_colormHG, *eb, VectorHG(1., 1., 1.)); // here white color
    // sample/nothing to do with curvature - create property map to store edge
    // color

    // ---
    const std::pair< BaseAdapterVisu *, HalfedgeGraph * > p =
        std::make_pair(_adapter, _mesh);

    if((*value_forceCompute) || (map_v_cmHG.find(p) == map_v_cmHG.end()))
    {
      auto m_gpm =
          get_property_map(FEVV::mesh_guiproperties, *_mesh, *pmaps_bag);

      map_v_cmHG[p] = std::make_tuple(
          FEVV::make_vertex_property_map< HalfedgeGraph,
                                          Filters::v_Curv< HalfedgeGraph > >(
              *_mesh),
          std::array< double, 4 >());

      curvature(_mesh,
                m_gpm,
                std::get< 0 >(map_v_cmHG[p]),
                std::get< 1 >(map_v_cmHG[p])[0],
                std::get< 1 >(map_v_cmHG[p])[1],
                std::get< 1 >(map_v_cmHG[p])[2],
                std::get< 1 >(map_v_cmHG[p])[3]);

#ifdef DBG_applyHG
      {
        VertexCurvatureMapHG curvMap = std::get< 0 >(map_v_cmHG[p]);
        uint v_count = 0;
        auto vit = vertices(*_mesh).first;
        auto vit_end = vertices(*_mesh).second;
        for(; vit != vit_end; ++vit)
        {
          v_count++;
          auto curvData = get(curvMap, *vit);
          std::cout << "vertex #" << v_count << std::endl;
          std::cout << "KmaxCurv = " << curvData.KmaxCurv << std::endl;
          std::cout << "KminCurv = " << curvData.KminCurv << std::endl;
          std::cout << "VKmaxCurv = " << curvData.VKmaxCurv[0] << " "
                    << curvData.VKmaxCurv[1] << " " << curvData.VKmaxCurv[2]
                    << std::endl;
          std::cout << "VKminCurv = " << curvData.VKminCurv[0] << " "
                    << curvData.VKminCurv[1] << " " << curvData.VKminCurv[2]
                    << std::endl;
        }
      }
#endif

      QMessageBox::information(0, "", QObject::tr("Curvature calculated."));
    }

    constructColorMap(*_mesh,
                      std::get< 0 >(map_v_cmHG[p]),
                      v_colormHG,
                      std::get< 1 >(map_v_cmHG[p])[0],
                      std::get< 1 >(map_v_cmHG[p])[1],
                      std::get< 1 >(map_v_cmHG[p])[2],
                      std::get< 1 >(map_v_cmHG[p])[3],
                      (*value_colorField)); // 1 : min, 2 : max
    // ---

    auto viewer = dynamic_cast< SimpleViewer * >(_adapter->getViewer());
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
#endif

#ifdef FEVV_USE_AIF
  void apply(BaseAdapterVisu *_adapter,
             MeshAIF *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    applyHG< MeshAIF >(_adapter, _mesh, pmaps_bag);
  }
#endif


  QStringList Generic_plugins() const override
  {
    return QStringList() << "CurvaturePlugin";
  }

  bool Generic_plugin(const QString &plugin) override
  {
    SimpleWindow *sw = static_cast< SimpleWindow * >(window);
      // dynamic_cast fails under OS X
    sw->onModificationParam("curvature_qt_p", this);
    sw->onApplyButton();

    return true;
  }

signals:
  void resetSignal();

protected:
  bool *value_forceCompute = new bool(false);

  bool *value_isGeod = new bool(false);
  double *value_radius = new double(0.0);

  int *value_colorField = new int(0);

  bool *value_displayMinDirections = new bool(false);
  bool *value_displayMaxDirections = new bool(false);

  // ---

  /*! \brief Look-up table for color rendering (from blue to red)*/
  double LUT_CourbureClust[3 * 256];
};

} // namespace FEVV

