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
#include "Dialogs/DialogDecompressionValence1.h"

#ifndef Q_MOC_RUN // MT : very important to avoid the error : ' Parse error at
                  // "BOOST_JOIN" ' -> (qt4 pb with boost)
#include "Visualization/PluginFilters/BasePluginQt.h"
#include "Visualization/SimpleViewer.h"

#include "Visualization/SimpleWindow.h"

#include "FEVV/Filters/Generic/Manifold/Compression_Valence/decompression_valence.h" // A) include the header of the filter corresponding to your operation

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

class DecompressionValencePlugin : public QObject,
                                   public Generic_PluginInterface,
                                   public BasePluginQt
{
  Q_OBJECT
  Q_INTERFACES(FEVV::Generic_PluginInterface)
#if(FEVV_USE_QT5) // see at the end of .cpp for QT4
  Q_PLUGIN_METADATA(IID "DecompressionValencePlugin")
#endif

  /*public:
      using BasePlugin::apply;*/
public:
  DecompressionValencePlugin() = default;
  ~DecompressionValencePlugin() = default;

public:
  void init() override { init(true, "example.p3d", -1, true, false); }

  void init(bool _forceCompute,
            const std::string &_p3dFilePath,
            int _stop_level,
            bool _write_info,
            bool _write_intermediate_meshes)
  {
    *value_forceCompute = _forceCompute;

    p3dFilePath = _p3dFilePath;
    stop_level = _stop_level;
    write_info = _write_info;
    write_intermediate_meshes = _write_intermediate_meshes;
    keep_intermediate_meshes = true;
    intermediate_meshes_void = nullptr;
    intermediate_vertexColorMaps_void = nullptr;
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
  void process(HalfedgeGraph *&_mesh, FEVV::PMapsContainer *&pmaps_bag)
  {
    std::cout << "Asking to apply DecompressionValence filter on file '"
              << p3dFilePath << "'!" << std::endl;

#if 0
        //TODO-elo DBG
        std::cout << "Decompression Valence Plugin DBG infos" << std::endl;
        std::cout << " * p3dFilePath=" << p3dFilePath << std::endl;
        std::cout << " * write_info=" << write_info << std::endl;
        std::cout << " * write_intermediate_meshes=" << write_intermediate_meshes << std::endl;
        std::cout << " * keep_intermediate_meshes=" << keep_intermediate_meshes << std::endl;
#endif

    // retrieve geometry property map
    auto pm = get(boost::vertex_point, *_mesh);

    // create vertex color property map
    using VertexColorMap =
        typename FEVV::PMap_traits< FEVV::vertex_color_t,
                                    HalfedgeGraph >::pmap_type;
    VertexColorMap v_cm;

    // create intermediate meshes storage
    std::vector< HalfedgeGraph * > *intermediate_meshes = nullptr;
    std::vector< VertexColorMap * > *intermediate_vertexColorMaps = nullptr;
    if(keep_intermediate_meshes)
    {
      intermediate_meshes = new std::vector< HalfedgeGraph * >;
      intermediate_vertexColorMaps = new std::vector< VertexColorMap * >;
    }

    // apply Decompression Valence filter
    std::string result;
    try
    {
      result =
          FEVV::Filters::decompression_valence(*_mesh,
                                               &pm,
                                               &v_cm,
                                               p3dFilePath,
                                               write_info,
                                               intermediate_meshes,
                                               intermediate_vertexColorMaps,
                                               stop_level,
                                               write_intermediate_meshes);

      // existing property maps are no more valid due to topological changes
      // purge the property maps bag except the vertex color map
      FEVV::PMapsContainer new_pmaps_bag;
      FEVV::put_property_map(FEVV::vertex_color, *_mesh, new_pmaps_bag, v_cm);
      *pmaps_bag = new_pmaps_bag;
    }
    catch(std::runtime_error &e)
    {
      std::cout << e.what() << std::endl;
    }

    intermediate_meshes_void = static_cast< void * >(intermediate_meshes);
    intermediate_vertexColorMaps_void =
        static_cast< void * >(intermediate_vertexColorMaps);

    std::string message("Decompression completed.\n\n");
    message.append(result);
    QMessageBox::information(0, "", QString::fromStdString(message));
  }


  template< typename HalfedgeGraph >
  void applyHG(BaseAdapterVisu *_adapter)
  {
    // get filter parameters from dialog window
    DialogDecompressionValence1 dial1;
    dial1.setDecompressionValenceParams(p3dFilePath,
                                        stop_level,
                                        write_info,
                                        write_intermediate_meshes,
                                        keep_intermediate_meshes);
    if(dial1.exec() == QDialog::Accepted)
      dial1.getDecompressionValenceParams(p3dFilePath,
                                          stop_level,
                                          write_info,
                                          write_intermediate_meshes,
                                          keep_intermediate_meshes);
    else
      return; // abort applying filter

    // create a new mesh and its property-maps bag for the filter
    auto mesh = new HalfedgeGraph;
    auto pmaps_bag = new FEVV::PMapsContainer;

    // apply filter
    process(mesh, pmaps_bag);

    // draw mesh(es)
    SimpleViewer *viewer =
        dynamic_cast< SimpleViewer * >(_adapter->getViewer());
    if(viewer)
    {
      using VertexColorMap =
          typename FEVV::PMap_traits< FEVV::vertex_color_t,
                                      HalfedgeGraph >::pmap_type;

      // draw the uncompressed mesh
      viewer->draw_or_redraw_mesh(mesh, pmaps_bag, false, false, "uncompressed");

      // space_time mode ON
      viewer->m_space_time = true;

      // draw all intermediate meshes
      auto intermediate_meshes =
          static_cast< std::vector< HalfedgeGraph * > * >(
              intermediate_meshes_void);
      auto intermediate_vertexColorMaps =
          static_cast< std::vector< VertexColorMap * > * >(
              intermediate_vertexColorMaps_void);

      if(intermediate_meshes_void && intermediate_vertexColorMaps_void)
      {
        // loop over intermediate meshes
        for(int i = intermediate_meshes->size() - 1; i >=0; i--)
        {
          HalfedgeGraph *mesh_i = (*intermediate_meshes)[i];
          VertexColorMap *v_cm_i = (*intermediate_vertexColorMaps)[i];

          FEVV::PMapsContainer *pmaps_bag_i = new FEVV::PMapsContainer;
          FEVV::put_property_map(
              FEVV::vertex_color, *mesh_i, *pmaps_bag_i, *v_cm_i);

          // draw intermediate mesh
          viewer->draw_or_redraw_mesh(mesh_i,
                                      pmaps_bag_i,
                                      false,
                                      false,
                                      std::string("level") + std::to_string(i));
        }
      }

      delete(intermediate_meshes);
      intermediate_meshes_void = nullptr;
      delete(intermediate_vertexColorMaps);
      intermediate_vertexColorMaps_void = nullptr;
    }

    // TODO-elo-to-keep-parameters-between-calls   reset();

    viewer->activate_time_mode();

    viewer->frame();
  }


#ifdef FEVV_USE_CGAL
  void apply(BaseAdapterVisu *_adapter,
             MeshLCC *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    apply(_adapter, static_cast< void * >(nullptr), nullptr);
  }

  void apply(BaseAdapterVisu *_adapter,
             MeshSurface *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    apply(_adapter, static_cast< void * >(nullptr), nullptr);
  }

  void apply(BaseAdapterVisu *_adapter,
             MeshPolyhedron *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    apply(_adapter, static_cast< void * >(nullptr), nullptr);
  }
#endif


#ifdef FEVV_USE_OPENMESH
  void apply(BaseAdapterVisu *_adapter,
             MeshOpenMesh *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    apply(_adapter, static_cast< void * >(nullptr), nullptr);
  }
#endif


#ifdef FEVV_USE_AIF
  void apply(BaseAdapterVisu *_adapter,
             MeshAIF *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    apply(_adapter, static_cast< void * >(nullptr), nullptr);
  }
#endif


  // case where the plugin is applied when no mesh is opened
  void apply(BaseAdapterVisu *_adapter,
             void *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    // ask the user for the datastructure
    std::string mesh_type = chooseDatastructureMsgBox();
    if(mesh_type == "NONE")
      return; // cancel pressed, aborting

    // apply plugin
#ifdef FEVV_USE_CGAL
    if(mesh_type == "POLYHEDRON")
    {
      applyHG< MeshPolyhedron >(_adapter);
    }
    else if(mesh_type == "SURFACEMESH")
    {
      applyHG< MeshSurface >(_adapter);
    }
    else if(mesh_type == "LCC")
    {
      applyHG< MeshLCC >(_adapter);
    }
    else if(mesh_type == "CGALPOINTSET")
    {
      QMessageBox::information(
          0,
          "",
          QObject::tr(
              "Decompression Valence filter is not yet compatible with CGALPointSet!"));
    }
#endif

#ifdef FEVV_USE_OPENMESH
    if(mesh_type == "OPENMESH")
    {
      applyHG< MeshOpenMesh >(_adapter);
    }
#endif

#ifdef FEVV_USE_AIF
    if(mesh_type == "AIF")
    {
      QMessageBox::information(
          0,
          "",
          QObject::tr(
              "Decompression Valence filter is not yet compatible with AIF!"));
    }
#endif

#ifdef FEVV_USE_PCL
    if(mesh_type == "PCLPOINTCLOUD")
    {
      QMessageBox::information(
          0,
          "",
          QObject::tr(
              "Decompression Valence filter is not yet compatible with PCLPointCloud!"));
    }
#endif
  }


  QStringList Generic_plugins() const override
  {
    return QStringList() << "DecompressionValencePlugin";
  }

  bool Generic_plugin(const QString &plugin) override
  {
    SimpleWindow *sw = static_cast< SimpleWindow * >(window);
      // dynamic_cast fails under OS X
    sw->onModificationParam("decompressionvalence_qt_p", this);
    sw->onApplyButton();

    return true;
  }

signals:
  void resetSignal();

protected:
  bool *value_forceCompute = new bool(false);

  // filter parameters
  std::string p3dFilePath;
  int stop_level;
  bool write_info;
  bool keep_intermediate_meshes;  // keep into RAM
  bool write_intermediate_meshes; // write to files
  void *intermediate_meshes_void;
  void *intermediate_vertexColorMaps_void;
};

} // namespace FEVV

