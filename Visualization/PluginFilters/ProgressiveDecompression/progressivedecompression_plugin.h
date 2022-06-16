// Copyright (c) 2012-2022 University of Lyon and CNRS (France).
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
#include "Dialogs/progressivedecompression_dialog.h"

#ifndef Q_MOC_RUN // MT : very important to avoid the error : ' Parse error at
                  // "BOOST_JOIN" ' -> (qt4 pb with boost)
#include "Visualization/PluginFilters/BasePluginQt.h"
#include "Visualization/SimpleViewer.h"

#include "Visualization/SimpleWindow.h"

// A) include the header of the filter corresponding to your operation
#include "FEVV/Filters/CGAL/Progressive_Compression/progressive_decompression_filter.hpp"
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

class ProgressiveDecompressionPlugin : public QObject,
                         public Generic_PluginInterface,
                         public BasePluginQt
{
  Q_OBJECT
  Q_INTERFACES(FEVV::Generic_PluginInterface)
#if(FEVV_USE_QT5) // see at the end of .cpp for QT4
  Q_PLUGIN_METADATA(IID "ProgressiveDecompressionPlugin")
#endif

  /*public:
      using BasePlugin::apply;*/
public:
  ProgressiveDecompressionPlugin() = default;
  ~ProgressiveDecompressionPlugin() = default;

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
  void process(HalfedgeGraph& _mesh, FEVV::PMapsContainer& pmaps_bag)
  {
    std::cout << "Asking to ProgressiveDecompression mesh ! " << std::endl;
    
    // retrieve or create vertex-color property map
    using VertexColorMap =
        typename FEVV::PMap_traits< FEVV::vertex_color_t,
		                            HalfedgeGraph >::pmap_type;
    VertexColorMap v_cm;
    if(has_map(pmaps_bag, FEVV::vertex_color))
    {
      std::cout << "use existing vertex-color map" << std::endl;
      v_cm = get_property_map(FEVV::vertex_color, _mesh, pmaps_bag);
    }
    else
    {
      std::cout << "create vertex-color map" << std::endl;
      v_cm = make_property_map(FEVV::vertex_color, _mesh);
      // store property map in property maps bag
      put_property_map(FEVV::vertex_color, _mesh, pmaps_bag, v_cm);
    }
  
    // retrieve or create vertex-normal property map
    using VertexNormalMap =
        typename FEVV::PMap_traits< FEVV::vertex_normal_t,
		                            HalfedgeGraph >::pmap_type;
    VertexNormalMap v_nm;
    if(has_map(pmaps_bag, FEVV::vertex_normal))
    {
      std::cout << "use existing vertex-normal map" << std::endl;
      v_nm = get_property_map(FEVV::vertex_normal, _mesh, pmaps_bag);
    }
    else
    {
      std::cout << "create vertex-normal map" << std::endl;
      v_nm = make_property_map(FEVV::vertex_normal, _mesh);
      // store property map in property maps bag
      put_property_map(FEVV::vertex_normal, _mesh, pmaps_bag, v_nm);
    }


	typename FEVV::PMap_traits< FEVV::edge_color_t, HalfedgeGraph >::pmap_type e_cm;
    if(has_map(pmaps_bag, FEVV::edge_color))
    {
      std::cout << "using edge color property map" << std::endl;
      e_cm = get_property_map(FEVV::edge_color, _mesh, pmaps_bag);
      // e_cm_ptr = &e_cm;
    }
    else
    {
      std::cout << "create vertex-color map" << std::endl;
      e_cm = make_property_map(FEVV::edge_color, _mesh);
      // store property map in property maps bag
      put_property_map(FEVV::edge_color, _mesh, pmaps_bag, e_cm);
    }

    // retrieve point property map (aka geometry)
    auto pm = get(boost::vertex_point, _mesh);

    // apply filter
    // B) call the filter corresponding to your operation
    std::string result;
    try {
      FEVV::Filters::progressive_decompression_filter(_mesh, 
                                                      pm, 
                                                      v_cm, 
                                                      //v_nm, 
                                                      e_cm, 
                                                      mesh_binary_path, 
                                                      true);
    }
    catch (std::exception& e) {
      result = e.what();
    }
	
    // existing property maps are no more valid due to topological changes
    // purge the property maps bag
    FEVV::PMapsContainer clean_pmaps_bag;
    pmaps_bag = clean_pmaps_bag;

    std::string message("Progressive decompression completed.\n\n");
    message.append(result);
    QMessageBox::information(0, "", QString::fromStdString(message));
  }

  template< typename HalfedgeGraph >
  void applyHG(BaseAdapterVisu *_adapter)
  {
    // get filter parameters from dialog window
    ProgressiveDecompressionDialog dialog;
    //dialog.setParameters(value_x, value_y, value_z);
    if(dialog.exec() == QDialog::Accepted)
      dialog.getParameters(mesh_binary_path);
    else
      return; // abort applying filter
  
  
    // create a new mesh and its property-maps bag for the filter
    HalfedgeGraph* mesh = new HalfedgeGraph();
    FEVV::PMapsContainer* pmaps_bag = new FEVV::PMapsContainer();
  
    // apply filter
    process(*mesh, *pmaps_bag);

    SimpleViewer *viewer =
        dynamic_cast< SimpleViewer * >(_adapter->getViewer());
    if(viewer)
      viewer->draw_or_redraw_mesh(mesh, pmaps_bag, false, false, "uncompressed");

    // comment next line to keep parameters values between calls
    reset();

#if(FEVV_USE_QT5)
    // empty
#else
    viewer->frame(); // necessary or not ?
#endif
  }

#if FEVV_USE_OPENMESH
  void apply(BaseAdapterVisu* _adapter, MeshOpenMesh* /*_mesh*/, FEVV::PMapsContainer * /*pmaps_bag*/) override
  {
    apply(_adapter, static_cast<void *>(nullptr), nullptr);
  }
#endif

#ifdef FEVV_USE_CGAL
#if 1
  void apply(BaseAdapterVisu* _adapter, MeshLCC* /*_mesh*/, FEVV::PMapsContainer * /*pmaps_bag*/) override
  {
    apply(_adapter, static_cast<void *>(nullptr), nullptr);
  }
#endif
  void apply(BaseAdapterVisu* _adapter, MeshSurface* /*_mesh*/, FEVV::PMapsContainer * /*pmaps_bag*/) override
  {
    apply(_adapter, static_cast<void *>(nullptr), nullptr);
  }

  void apply(BaseAdapterVisu* _adapter, MeshPolyhedron* /*_mesh*/, FEVV::PMapsContainer * /*pmaps_bag*/) override
  {
    apply(_adapter, static_cast<void *>(nullptr), nullptr);
  }
#endif

#if FEVV_USE_AIF
#if 0 // halfedge iteration is not compatible with AIF yet (halfedge_iterator type cannot be exported)
  void apply(BaseAdapterVisu* _adapter, MeshAIF* /*_mesh*/, FEVV::PMapsContainer * /*pmaps_bag*/) override
  {
    apply(_adapter, static_cast<void *>(nullptr), nullptr);
  }
#endif
#endif

  // case where the plugin is applied when no mesh is opened
  void apply(BaseAdapterVisu *_adapter,
              void * /*_mesh*/,
              FEVV::PMapsContainer * /*pmaps_bag*/) override
  {
    // ask the user for the datastructure
    std::string mesh_type = chooseDatastructureMsgBox();
    if (mesh_type == "NONE")
      return; // cancel pressed, aborting

    // apply plugin
#ifdef FEVV_USE_CGAL
    if (mesh_type == "POLYHEDRON")
    {
      applyHG< MeshPolyhedron >(_adapter);
    }
    else if (mesh_type == "SURFACEMESH")
    {
      applyHG< MeshSurface >(_adapter);
    }
    else if (mesh_type == "LCC")
    {
      applyHG< MeshLCC >(_adapter);
    }
#endif

#ifdef FEVV_USE_OPENMESH
    if (mesh_type == "OPENMESH")
    {
		applyHG< MeshOpenMesh >(_adapter);
    }
#endif

#ifdef FEVV_USE_AIF
    if (mesh_type == "AIF")
    {
      QMessageBox::information(
        0,
        "",
        QObject::tr(
          "Progressive Decompression filter is not yet compatible with AIF!"));
    }
#endif
  }

  QStringList Generic_plugins() const override
  {
    return QStringList() << "ProgressiveDecompressionPlugin";
  }

  bool Generic_plugin(const QString &/*plugin*/) override
  {
    SimpleWindow *sw = static_cast< SimpleWindow * >(window);
      // dynamic_cast fails under OS X
    sw->onModificationParam("progressivedecompression_qt_p", this);
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
  std::string mesh_binary_path ="";
};

} // namespace FEVV

