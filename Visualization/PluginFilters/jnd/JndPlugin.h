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
#include "Dialogs/DialogJnd1.h"

#ifndef Q_MOC_RUN // MT : very important to avoid the error : ' Parse error at
                  // "BOOST_JOIN" ' -> (qt4 pb with boost)
#include <sstream>
#include <string>
#include <boost/foreach.hpp>
#include <time.h>
#include <cmath>

#include "Visualization/PluginFilters/BasePluginQt.h"
#include "Visualization/SimpleViewer.h"

#include "Visualization/SimpleWindow.h"

#include "FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/just_noticeable_distortion.hpp" // A) include the header of the filter corresponding to your operation


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

class JndPlugin : public QObject,
                  public Generic_PluginInterface,
                  public BasePluginQt
{
  Q_OBJECT
  Q_INTERFACES(FEVV::Generic_PluginInterface)
#if(FEVV_USE_QT5) // see at the end of .cpp for QT4
  Q_PLUGIN_METADATA(IID "JndPlugin")
#endif

  /*public:
      using BasePlugin::apply;*/
public:
  JndPlugin() = default;
  ~JndPlugin() = default;

public:
  void init() override
  {
    *screen_width = 1920;
    *screen_height = 1080;
    *screen_size = 55.;
    *user_dist = 50.;
    *scene_height = 1080;
    *scene_fov = M_PI * 0.3333;
    *number_of_lights = 128;
    *use_log = true;
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
  void process(HalfedgeGraph *_m, FEVV::PMapsContainer *pmaps_bag)
  {

    typedef boost::graph_traits< HalfedgeGraph > GraphTraits;
    typedef typename GraphTraits::vertex_iterator vertex_iterator;
    typedef typename GraphTraits::vertex_descriptor vertex_descriptor;

    ScreenParam screen(*screen_width, *screen_height, *screen_size);
    UserParam user(*user_dist);
    SceneParam scene(*scene_height, *scene_fov);


    // Note: the property maps must be extracted from the
    //       property maps bag, and explicitely passed as
    //       parameters to the filter, in order to make
    //       clear what property is used by the filter

    // retrieve point property map (aka geometry)
    auto pm = get(boost::vertex_point, *_m);

    // retrieve or create vertex-color property map
    using VertexColorMap =
        typename FEVV::PMap_traits< FEVV::vertex_color_t,
                                    HalfedgeGraph >::pmap_type;
    typedef typename boost::property_traits< VertexColorMap >::value_type Color;
    VertexColorMap v_cm;
    if(has_map(*pmaps_bag, FEVV::vertex_color))
    {
      std::cout << "use existing vertex-color map" << std::endl;
      v_cm = get_property_map(FEVV::vertex_color, *_m, *pmaps_bag);
    }
    else
    {
      std::cout << "create vertex-color map" << std::endl;
      v_cm = make_property_map(FEVV::vertex_color, *_m);
      // store property map in property maps bag
      put_property_map(FEVV::vertex_color, *_m, *pmaps_bag, v_cm);
    }

    // retrieve or create vertex-normal property map
    using FaceNormalMap =
        typename FEVV::PMap_traits< FEVV::face_normal_t,
                                    HalfedgeGraph >::pmap_type;
    FaceNormalMap f_nm;
    if(has_map(*pmaps_bag, FEVV::face_normal))
    {
      std::cout << "use existing face-normal map" << std::endl;
      f_nm = get_property_map(FEVV::face_normal, *_m, *pmaps_bag);
    }
    else
    {
      std::cout << "create face-normal map" << std::endl;
      f_nm = make_property_map(FEVV::face_normal, *_m);
      // store property map in property maps bag
      put_property_map(FEVV::face_normal, *_m, *pmaps_bag, f_nm);
      FEVV::Filters::calculate_face_normals(*_m, pm, f_nm);
    }

    // retrieve or create vertex-normal property map
    using VertexNormalMap =
        typename FEVV::PMap_traits< FEVV::vertex_normal_t,
                                    HalfedgeGraph >::pmap_type;
    VertexNormalMap v_nm;
    if(has_map(*pmaps_bag, std::string("v:normal")))
    {
      std::cout << "use existing vertex-normal map" << std::endl;
      v_nm = boost::any_cast< VertexNormalMap >(pmaps_bag->at("v:normal"));
    }
    else
    {
      std::cout << "create vertex-normal map" << std::endl;
      v_nm = make_property_map(FEVV::vertex_normal, *_m);
      // store property map in property maps bag
      put_property_map(FEVV::vertex_normal, *_m, *pmaps_bag, v_nm);
      FEVV::Filters::calculate_vertex_normals(*_m, pm, f_nm, v_nm);
    }

    using JndTypeMap =
        typename FEVV::Vertex_pmap_traits< HalfedgeGraph, double >::pmap_type;
    JndTypeMap jnd_m;
    clock_t tStart = clock();

    if(has_map(*pmaps_bag, std::string("v:jnd")) && !*force_jnd)
    {
      std::cout << "use existing jnd map" << std::endl;
      jnd_m = boost::any_cast< JndTypeMap >(pmaps_bag->at("v:jnd"));
    }
    else
    {
      jnd_m = FEVV::make_vertex_property_map< HalfedgeGraph, double >(*_m);
      std::cout << "computing jnd map" << std::endl;
      FEVV::Filters::just_noticeable_distortion_filter(
          *_m, pm, v_nm, f_nm, jnd_m, screen, user, scene, *number_of_lights);
      (*pmaps_bag)["v:jnd"] = jnd_m;
    }
    clock_t tJND = clock();

    double min_jnd = 10000.0;
    double max_jnd = 0.0;

    BOOST_FOREACH(vertex_descriptor vi, vertices(*_m))
    {
      auto jnd = get(jnd_m, vi);
      if(min_jnd > jnd)
        min_jnd = jnd;
      if(max_jnd < jnd)
        max_jnd = jnd;
    }
    clock_t tMM = clock();

    BOOST_FOREACH(vertex_descriptor vi, vertices(*_m))
    {
      auto jnd = get(jnd_m, vi);
      double color;
      if(*use_log)
        color = 2.0 * ((std::log((*log_disp + jnd))) /
                       (std::log((*log_disp + max_jnd)))) -
                1.0;
      else
        color = 2.0 * ((jnd - min_jnd) / (max_jnd - min_jnd)) - 1.0;
      Color newcolor(red(color), green(color), blue(color));
      put(v_cm, vi, newcolor);
    }

    std::cout << "Time used to compute JND : "
              << (double)(tJND - tStart) / CLOCKS_PER_SEC << std::endl;
    std::cout << "Time used to compute minmax : "
              << (double)(tMM - tJND) / CLOCKS_PER_SEC << std::endl;
    std::cout << "Time used to compute colors : "
              << (double)(clock() - tMM) / CLOCKS_PER_SEC << std::endl;
    std::cout << "JND min : " << min_jnd << std::endl;
    std::cout << "JND max : " << max_jnd << std::endl;
  }


  template< typename HalfedgeGraph >
  void applyHG(BaseAdapterVisu *_adapter,
               HalfedgeGraph *_mesh,
               FEVV::PMapsContainer *pmaps_bag)
  {
    // get filter parameters from dialog window
    DialogJnd1 dial1;
    if(dial1.exec() == QDialog::Accepted)
    {
      dial1.getProcess(*screen_width,
                       *screen_height,
                       *screen_size,
                       *scene_height,
                       *scene_fov,
                       *user_dist,
                       *number_of_lights,
                       *log_disp,
                       *use_log,
                       *force_jnd);

      *scene_fov = M_PI * *scene_fov;
    }
    else
      return; // abort applying filter

    // apply filter
    process(_mesh, pmaps_bag);

    // redraw mesh
    SimpleViewer *viewer =
        dynamic_cast< SimpleViewer * >(_adapter->getViewer());

    if(viewer)
    {
      viewer->draw_or_redraw_mesh(_mesh, pmaps_bag, true, false);
    }

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
    return QStringList() << "JndPlugin";
  }

  bool Generic_plugin(const QString &plugin) override
  {
    SimpleWindow *sw = static_cast< SimpleWindow * >(window);
      // dynamic_cast fails under OS X
    sw->onModificationParam("jnd_qt_p", this);
    sw->onApplyButton();

    return true;
  }

signals:
  void resetSignal();

protected:
  int *screen_width = new int(1920);
  int *screen_height = new int(1080);
  double *screen_size = new double(55.);
  double *user_dist = new double(50.);
  int *scene_height = new int(1080);
  double *scene_fov = new double(M_PI * 0.3333);
  int *number_of_lights = new int(128);
  float *log_disp = new float(0.5);
  bool *use_log = new bool(true);
  bool *force_jnd = new bool(true);

  double interpolate(double val, double y0, double x0, double y1, double x1)
  {
    return (val - x0) * (y1 - y0) / (x1 - x0) + y0;
  }

  double base(double val)
  {
    if(val >= 1.0)
      return 1.0;
    else if(val > 0.0)
      return (1.0 - val);
    else
      return 0.0;
  }

  double red(double gray) { return base(1.0 + gray); }
  double green(double gray) { return base(1.0 - fabs(gray)); }
  double blue(double gray) { return base(1.0 - gray); }
};

} // namespace FEVV

