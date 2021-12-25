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
#include "Dialogs/progressivecompression_dialog.h"

#ifndef Q_MOC_RUN // MT : very important to avoid the error : ' Parse error at
                  // "BOOST_JOIN" ' -> (qt4 pb with boost)
#include "Visualization/PluginFilters/BasePluginQt.h"
#include "Visualization/SimpleViewer.h"

#include "Visualization/SimpleWindow.h"

// A) include the header of the filter corresponding to your operation
#include "FEVV/Filters/CGAL/Progressive_Compression/progressive_compression_filter.hpp"
#include "FEVV/Wrappings/properties.h"

#ifdef FEVV_USE_CGAL
#include "FEVV/Wrappings/properties_polyhedron_3.h"
#include "FEVV/Wrappings/properties_surface_mesh.h"
#include "FEVV/Wrappings/properties_linear_cell_complex.h"
#endif // FEVV_USE_CGAL
#if FEVV_USE_OPENMESH
#include "FEVV/Wrappings/properties_openmesh.h"
#endif // FEVV_USE_OPENMESH
#if FEVV_USE_AIF
#include "FEVV/Wrappings/properties_aif.h"
#endif // FEVV_USE_AIF
#endif

namespace FEVV {

	template< typename MeshT >
void
set_mesh_and_properties(
    MeshT &m,
    FEVV::PMapsContainer &pmaps_bag,
    typename FEVV::PMap_traits< FEVV::vertex_color_t, MeshT >::pmap_type &v_cm,
    typename FEVV::PMap_traits< FEVV::edge_color_t, MeshT >::pmap_type &e_cm,
    typename FEVV::PMap_traits< FEVV::face_normal_t, MeshT >::pmap_type &/*f_nm*/,

    typename FEVV::PMap_traits< FEVV::vertex_normal_t, MeshT >::pmap_type &v_nm)

{
  // Note: the property maps must be extracted from the
  //       property maps bag, and explicitely passed as
  //       parameters to the filter, in order to make
  //       clear what property is used by the filter

  // retrieve or create vertex-color property map
  if(has_map(pmaps_bag, FEVV::vertex_color))
  {
    std::cout << "use existing vertex-color map" << std::endl;
    v_cm = get_property_map(FEVV::vertex_color, m, pmaps_bag);
  }
  else
  {
    std::cout << "create vertex-color map" << std::endl;
    v_cm = make_property_map(FEVV::vertex_color, m);
    // store property map in property maps bag
    put_property_map(FEVV::vertex_color, m, pmaps_bag, v_cm);
  }
  if(has_map(pmaps_bag, FEVV::edge_color))
  {
    std::cout << "using edge color property map" << std::endl;
    e_cm = get_property_map(FEVV::edge_color, m, pmaps_bag);
    // e_cm_ptr = &e_cm;
  }
  else
  {
    std::cout << "create vertex-color map" << std::endl;
    e_cm = make_property_map(FEVV::edge_color, m);
    // store property map in property maps bag
    put_property_map(FEVV::edge_color, m, pmaps_bag, e_cm);
  }
  // retrieve or create vertex-normal property map

  if(has_map(pmaps_bag, FEVV::vertex_normal))
  {
    std::cout << "use existing vertex-normal map" << std::endl;
    v_nm = get_property_map(FEVV::vertex_normal, m, pmaps_bag);
  }
  else
  {
    std::cout << "create vertex-normal map" << std::endl;
    v_nm = make_property_map(FEVV::vertex_normal, m);
    // store property map in property maps bag
    put_property_map(FEVV::vertex_normal, m, pmaps_bag, v_nm);
  }
}
class ProgressiveCompressionPlugin : public QObject,
                         public Generic_PluginInterface,
                         public BasePluginQt
{
  Q_OBJECT
  Q_INTERFACES(FEVV::Generic_PluginInterface)
#if(FEVV_USE_QT5) // see at the end of .cpp for QT4
  Q_PLUGIN_METADATA(IID "ProgressiveCompressionPlugin")
#endif

  /*public:
      using BasePlugin::apply;*/
public:
  ProgressiveCompressionPlugin() = default;
  ~ProgressiveCompressionPlugin() = default;

public:
  void init() override {}

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
  void process(HalfedgeGraph *_mesh, FEVV::PMapsContainer *pmaps_bag)
  {
    std::cout << "Asking to ProgressiveCompression mesh ! " << std::endl;

    // retrieve or create vertex-color property map
    using VertexColorMap =
        typename FEVV::PMap_traits< FEVV::vertex_color_t,
		                            HalfedgeGraph >::pmap_type;
    VertexColorMap v_cm;
    
    // retrieve or create vertex-normal property map
    using VertexNormalMap =
        typename FEVV::PMap_traits< FEVV::vertex_normal_t,
		                            HalfedgeGraph >::pmap_type;
    VertexNormalMap v_nm;
    
	typename FEVV::PMap_traits<FEVV::edge_color_t, HalfedgeGraph>::pmap_type e_cm;
    typename FEVV::PMap_traits< FEVV::face_normal_t, HalfedgeGraph >::pmap_type f_nm;
    set_mesh_and_properties(*_mesh,
                            *pmaps_bag,

                            v_cm,
                            e_cm,
                            f_nm,
                            v_nm);
        

	typename FEVV::Geometry_traits< HalfedgeGraph >  gt(*_mesh); // retrieve point property map (aka
                                            // geometry)
    auto pm = get(boost::vertex_point, *_mesh);

    // apply filter
    // B) call the filter corresponding to your operation
    FEVV::Filters::Parameters params(_predictor,_operator,_metric,true,false,_quantization);
    std::string filepath = "null";
	
    std::string result;
    try {
      progressive_compression_filter(*_mesh,
                                     pm,
                                     v_cm,
                                     e_cm,
                                     //v_nm,
                                     gt,
                                     params,
                                     filepath,
                                     _filepath,
                                     _nb_batches,
                                     _minimum_vertices,
                                     _batch_stop,
                                     true,
                                     true,
                                     false,
                                     "progressive_compression_original_mesh_after_preprocess_plugin.off");									 
    }
    catch (std::exception& e) {
      result = e.what();
    }

    // existing property maps are no more valid due to topological changes
    // purge the property maps bag
    FEVV::PMapsContainer clean_pmaps_bag;
    *pmaps_bag = clean_pmaps_bag;

    std::string message("Progressive compression completed.\n\n");
    message.append(result);
    QMessageBox::information(0, "", QString::fromStdString(message));
  }

  template< typename HalfedgeGraph >
  void applyHG(BaseAdapterVisu *_adapter,
               HalfedgeGraph *_mesh,
               FEVV::PMapsContainer *pmaps_bag)
  {
    // get filter parameters from dialog window
    ProgressiveCompressionDialog dialog;
    std::string pred = "";
    std::string metr = "";
    std::string vkept = "";
    std::string filepath = "";
    std::string batch_string = "";
    if(dialog.exec() == QDialog::Accepted)
      dialog.getParameters(_quantization, _nb_batches,_minimum_vertices, metr, vkept, pred,filepath, batch_string);
    else
      return; // abort applying filter

	_filepath = filepath;
    _predictor = FEVV::Filters::PREDICTION_TYPE::BUTTERFLY;
        _batch_stop = FEVV::Filters::BATCH_CONDITION::ALL_EDGES;

	if(pred.compare("Butterfly") == 0)
    {
      _predictor = FEVV::Filters::PREDICTION_TYPE::BUTTERFLY;
    }
    if(pred.compare("Delta") == 0)
    {
      _predictor =  FEVV::Filters::PREDICTION_TYPE::DELTA;
    }
    

	_metric = FEVV::Filters::METRIC_TYPE::QEM;
	 if(metr.compare("QEM") == 0)
    {
      _metric = FEVV::Filters::METRIC_TYPE::QEM;
    }
    if(metr.compare("Volume Preserving") == 0)
    {
      _metric =  FEVV::Filters::METRIC_TYPE::VOLUME_PRESERVING;
    }
    if(metr.compare("EdgeLength") == 0)
    {
      _metric = FEVV::Filters::METRIC_TYPE::EDGE_LENGTH;
    }

	if(batch_string.compare("All Edges") == 0)
    {
          _batch_stop = FEVV::Filters::BATCH_CONDITION::ALL_EDGES;
    }
    if(batch_string.compare("Reach Threshold") == 0)
    {
      _batch_stop = FEVV::Filters::BATCH_CONDITION::REACH_THRESHOLD;
    }

	_operator = FEVV::Filters::VKEPT_POSITION::MIDPOINT;
	if(vkept.compare("MidPoint") == 0)
    {
      _operator = FEVV::Filters::VKEPT_POSITION::MIDPOINT;
    }
    if(vkept.compare("Halfedge") == 0)
    {
      _operator = FEVV::Filters::VKEPT_POSITION::HALFEDGE;
    }
   

    // apply filter
    process(_mesh, pmaps_bag);

    SimpleViewer *viewer =
        dynamic_cast< SimpleViewer * >(_adapter->getViewer());
    if(viewer)
      viewer->draw_or_redraw_mesh(_mesh, pmaps_bag, true, false);

    // comment next line to keep parameters values between calls
    reset();

    viewer->frame();
  }

#ifdef FEVV_USE_OPENMESH
#if 0 // Works with CGAL AABB, no solution available for OpenMesh yet
  void apply(BaseAdapterVisu *_adapter,
             MeshOpenMesh *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    applyHG< MeshOpenMesh >(_adapter, _mesh, pmaps_bag);
  }
#endif
#endif

#ifdef FEVV_USE_CGAL
#if 1
  void apply(BaseAdapterVisu *_adapter,
             MeshLCC *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    applyHG< MeshLCC >(_adapter, _mesh, pmaps_bag);
  }
#endif
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
#if 0 // halfedge iteration is not compatible with AIF yet (halfedge_iterator type cannot be exported)
      // Works with CGAL AABB, no solution available for AIF yet
  void apply(BaseAdapterVisu *_adapter,
             MeshAIF *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    applyHG< MeshAIF >(_adapter, _mesh, pmaps_bag);
  }
#endif
#endif

  QStringList Generic_plugins() const override
  {
    return QStringList() << "ProgressiveCompressionPlugin";
  }

  bool Generic_plugin(const QString &/*plugin*/) override
  {
    SimpleWindow *sw = static_cast< SimpleWindow * >(window);
      // dynamic_cast fails under OS X
    sw->onModificationParam("progressivecompression_qt_p", this);
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
  Filters::PREDICTION_TYPE _predictor;
  Filters::METRIC_TYPE _metric;
  Filters::VKEPT_POSITION _operator;
  Filters::BATCH_CONDITION _batch_stop;
  std::string _filepath;
  int _nb_batches;
  int _quantization;
  int _minimum_vertices;
};

} // namespace FEVV

