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
#include "Dialogs/DialogCompressionValence1.h"

#ifndef Q_MOC_RUN // MT : very important to avoid the error : ' Parse error at
                  // "BOOST_JOIN" ' -> (qt4 pb with boost)
#include "Visualization/PluginFilters/BasePluginQt.h"
#include "Visualization/SimpleViewer.h"

#include "Visualization/SimpleWindow.h"

#include "FEVV/Filters/Generic/Manifold/Compression_Valence/compression_valence.h" // A) include the header of the filter corresponding to your operation

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

class CompressionValencePlugin : public QObject,
                                 public Generic_PluginInterface,
                                 public BasePluginQt
{
  Q_OBJECT
  Q_INTERFACES(FEVV::Generic_PluginInterface)
#if(FEVV_USE_QT5) // see at the end of .cpp for QT4
  Q_PLUGIN_METADATA(IID "CompressionValencePlugin")
#endif

  /*public:
      using BasePlugin::apply;*/
public:
  CompressionValencePlugin() = default;
  ~CompressionValencePlugin() = default;

public:
  void init() override { init(true, "compressed_mesh.p3d", false, 10, 100); }

  void init(bool _forceCompute,
            const std::string &_p3dFilePath,
            bool _with_adaptative_quantization,
            int _quantization_bits,
            int _max_vertices)
  {
    *value_forceCompute = _forceCompute;

    p3dFilePath = _p3dFilePath;
    with_adaptative_quantization = _with_adaptative_quantization;
    quantization_bits = _quantization_bits;
    max_vertices = _max_vertices;
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
  void process(HalfedgeGraph *_mesh, FEVV::PMapsContainer *pmaps_bag)
  {
    std::cout << "Asking to apply CompressionValence filter ! " << std::endl;

    // do not compress if p3d file name has not been set
    bool with_compression = !p3dFilePath.empty();

#if 0
        //TODO-elo DBG
        std::cout << "Compression Valence Plugin DBG infos" << std::endl;
        std::cout << " * p3dFilePath=" << p3dFilePath << std::endl;
        std::cout << " * with_compression=" << with_compression << std::endl;
        std::cout << " * with_adaptative_quantization=" << with_adaptative_quantization << std::endl;
        std::cout << " * quantization_bits=" << quantization_bits << std::endl;
        std::cout << " * max_vertices=" << max_vertices << std::endl;
#endif

    // retrieve geometry property map
    auto pm = get(boost::vertex_point, *_mesh);

    // retrieve vertex color property map
    using VertexColorMap =
        typename FEVV::PMap_traits< FEVV::vertex_color_t,
                                    HalfedgeGraph >::pmap_type;
    VertexColorMap v_cm;
    VertexColorMap *v_cm_ptr = nullptr;

    if(pmaps_bag && has_map(*pmaps_bag, FEVV::vertex_color))
    {
      std::cout << "Compression Valence filter: using vertex color property map"
                << std::endl;
      v_cm = FEVV::get_property_map(FEVV::vertex_color, *_mesh, *pmaps_bag);
      v_cm_ptr = &v_cm;
    }
    else
      std::cout << "Compression Valence filter: no vertex color property map "
                   "found, ignoring color"
                << std::endl;

    // copy the original mesh and its vertex-color property map (if any)
    // to be able to display the original mesh and the compressed
    // mesh at the same time in space mode

    // create new empty mesh
    HalfedgeGraph *mesh_copy = new HalfedgeGraph;
    mesh_copy_void = static_cast< void * >(mesh_copy);

    // create new vertex-color map if needed
    VertexColorMap *v_cm_copy_ptr = nullptr;
    if(v_cm_ptr != nullptr)
    {
      v_cm_copy_ptr = new VertexColorMap;
    }

    // put vertex-color map into a property maps bag
    pmaps_bag_copy = new FEVV::PMapsContainer;
    if(v_cm_copy_ptr != nullptr)
    {
      FEVV::put_property_map(
          FEVV::vertex_color, *mesh_copy, *pmaps_bag_copy, *v_cm_copy_ptr);
    }

    // copy mesh and vertex-color
    using PointMap =
        typename boost::property_map< HalfedgeGraph,
                                      boost::vertex_point_t >::type;
    Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >::
        copy_mesh(*_mesh, v_cm_ptr, *mesh_copy, v_cm_copy_ptr);

    // apply Compression Valence filter on original mesh
    std::string result;
    result = FEVV::Filters::compression_valence(*_mesh,
                                                &pm,
                                                v_cm_ptr,
                                                "",
                                                p3dFilePath,
                                                with_compression,
                                                with_adaptative_quantization,
                                                max_vertices,
                                                quantization_bits);

    // existing property maps are no more valid due to topological changes
    // purge the property maps bag except the vertex color map
    FEVV::PMapsContainer clean_pmaps_bag;
    if(v_cm_ptr != nullptr)
    {
      // put vertex color property map into property maps bag
      FEVV::put_property_map(FEVV::vertex_color, *_mesh, clean_pmaps_bag, v_cm);
    }
    *pmaps_bag = clean_pmaps_bag;

    std::string message("Compression completed.\n\n");
    message.append(result);
    QMessageBox::information(0, "", QString::fromStdString(message));
  }

  template< typename HalfedgeGraph >
  void applyHG(BaseAdapterVisu *_adapter,
               HalfedgeGraph *_mesh,
               FEVV::PMapsContainer *pmaps_bag)
  {
    // get filter parameters from dialog window
    DialogCompressionValence1 dial1;
    dial1.setCompressionValenceParams(p3dFilePath,
                                      with_adaptative_quantization,
                                      quantization_bits,
                                      max_vertices);
    if(dial1.exec() == QDialog::Accepted)
    {
      dial1.getCompressionValenceParams(p3dFilePath,
                                        with_adaptative_quantization,
                                        quantization_bits,
                                        max_vertices);

      // if the P3D file name doesn't end with the right extension, fix it
      if(p3dFilePath.size() < 4 ||
         p3dFilePath.substr(p3dFilePath.size() - 4) != ".p3d")
        p3dFilePath.append(".p3d");
    }
    else
      return; // abort applying filter

    // apply filter
    process(_mesh, pmaps_bag);

    // redraw/draw meshes
    SimpleViewer *viewer =
        dynamic_cast< SimpleViewer * >(_adapter->getViewer());
    if(viewer)
    {
      // redraw main mesh -> which is NOW compressed
      viewer->draw_or_redraw_mesh(_mesh, pmaps_bag, true, true, "compressed");

      // space_time mode ON
      viewer->m_space_time = true;

      // draw old original mesh
      viewer->draw_or_redraw_mesh(
          static_cast< HalfedgeGraph * >(mesh_copy_void),
          pmaps_bag_copy,
          false,
          false,
          "original_mesh_copied");
    }

    reset();

    viewer->activate_time_mode();

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
#if 0 //TODO-elo  restore when Compression Valence compiles with AIF
  void apply(BaseAdapterVisu *_adapter,
             MeshAIF *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
		
		applyHG<MeshAIF>(_adapter, _mesh, pmaps_bag);
  }
#endif
#endif


  QStringList Generic_plugins() const override
  {
    return QStringList() << "CompressionValencePlugin";
  }

  bool Generic_plugin(const QString &plugin) override
  {
    SimpleWindow *sw = static_cast< SimpleWindow * >(window);
      // dynamic_cast fails under OS X
    sw->onModificationParam("compressionvalence_qt_p", this);
    sw->onApplyButton();

    return true;
  }

signals:
  void resetSignal();

protected:
  bool *value_forceCompute = new bool(false);

  // filter parameters
  std::string p3dFilePath;
  bool with_adaptative_quantization;
  int quantization_bits;
  int max_vertices;

  void *mesh_copy_void;
  FEVV::PMapsContainer *pmaps_bag_copy;
};

} // namespace FEVV

