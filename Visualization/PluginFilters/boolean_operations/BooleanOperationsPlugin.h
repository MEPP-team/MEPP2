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
#include "Dialogs/DialogBooleanOperations1.h"

#ifndef Q_MOC_RUN // MT : very important to avoid the error : ' Parse error at
                  // "BOOST_JOIN" ' -> (qt4 pb with boost)
#include "Visualization/PluginFilters/BasePlugin.h"
#include "Visualization/SimpleViewer.h"

#include "Visualization/SimpleWindow.h"

#include "FEVV/Filters/CGAL/Polyhedron/Boolean_Operations/boolean_operations.hpp"

#include "FEVV/Wrappings/properties.h"

// Polyhedron is used internally
#include "FEVV/Wrappings/properties_polyhedron_3.h"

#ifdef FEVV_USE_CGAL
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

class BooleanOperationsPlugin : public QObject,
                         public Generic_PluginInterface,
                         public BasePlugin
{
  Q_OBJECT
  Q_INTERFACES(FEVV::Generic_PluginInterface)
#if(FEVV_USE_QT5) // see at the end of .cpp for QT4
  Q_PLUGIN_METADATA(IID "BooleanOperationsPlugin")
#endif

  /*public:
      using BasePlugin::apply;*/
public:
  BooleanOperationsPlugin() = default;
  ~BooleanOperationsPlugin() = default;

public:
  void init() override
  {
    m_operation = "UNION";
    m_output_mesh_void = nullptr;
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
  void process(HalfedgeGraph *mesh_A, HalfedgeGraph *mesh_B)
  {
    std::cout << "Asking to apply BooleanOperations filter ! " << std::endl;

    // create output mesh
    HalfedgeGraph *output_mesh = new HalfedgeGraph;

    // apply filter
    if(m_operation == "UNION")
      FEVV::Filters::boolean_union(*mesh_A, *mesh_B, *output_mesh);
    else if(m_operation == "INTER")
      FEVV::Filters::boolean_inter(*mesh_A, *mesh_B, *output_mesh);
    else
      FEVV::Filters::boolean_minus(*mesh_A, *mesh_B, *output_mesh);

    // store output mesh for later display
    m_output_mesh_void = static_cast< void * >(output_mesh);
  }

  template< typename HalfedgeGraph >
  void applyHG(BaseAdapterVisu *_adapter,
               HalfedgeGraph *_mesh,
               FEVV::PMapsContainer *pmaps_bag)
  {
    SimpleViewer< HalfedgeGraph > *viewer =
      dynamic_cast< SimpleViewer< HalfedgeGraph > * >(_adapter->getViewer());

    // retrieve the two input meshes in current viewer window,
    // then apply the filter

    std::vector< HalfedgeGraph * > meshes = viewer->getMeshes();
    if(meshes.size() >= 2)
    {
      auto mA = meshes[0];
      auto mB = meshes[1];

      process(mA, mB); // apply filter
    }
    else
    {
      QMessageBox::information(0,
          "",
          QObject::tr("Boolean Operations filter "
            "needs two meshes "
            "opened in Space or Time."));
    }

    // draw output mesh
    if(viewer)
    {
      // space_time mode ON
      viewer->m_space_time = true;

      // draw output mesh
      auto output_mesh = static_cast< HalfedgeGraph * >( m_output_mesh_void);
      if(output_mesh)
      {
          // pmaps_bag is required for display
          auto output_pmaps_bag = new FEVV::PMapsContainer;
          viewer->draw_or_redraw_mesh(output_mesh,
                                      output_pmaps_bag,
                                      false,
                                      false,
                                      m_operation);
      }
    }

    //ELO comment next line to keep parameters between calls
    //reset();

    viewer->frame();
  }

#ifdef FEVV_USE_OPENMESH
  void apply(BaseAdapterVisu *_adapter,
             MeshOpenMesh *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
#if 0
		//ELO-note: compiling error in Cartesian_converter.h with OM
    applyHG< MeshOpenMesh >(_adapter, _mesh, pmaps_bag);
#else
    QMessageBox::information(
        0,
        "",
        QObject::tr(
            "Boolean Operations filter is not yet compatible with OpenMesh!"));
#endif
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
#if 0
		//ELO-note: missing num_halfedges() used by CGAL::copy_face_graph()
    //          and compiling error in Cartesian_converter.h with AIF
		applyHG<MeshAIF>(_adapter, _mesh, pmaps_bag);
#else
    QMessageBox::information(
        0,
        "",
        QObject::tr(
            "Boolean Operations filter is not yet compatible with AIF!"));
#endif
  }
#endif

  QStringList Generic_plugins() const override
  {
    return QStringList() << "BooleanOperationsPlugin";
  }

  bool Generic_plugin(const QString &plugin) override
  {
    DialogBooleanOperations1 dial1;
    dial1.setParameters(m_operation);
    if(dial1.exec() == QDialog::Accepted)
    {
      dial1.getParameters(m_operation);

      SimpleWindow *sw = static_cast< SimpleWindow * >(window);
          // dynamic_cast fails under OS X

      sw->onModificationParam("booleanoperations_qt_p", this);
      sw->onApplyButton();

      return true;
    }

    return false;
  }

signals:
  void resetSignal();

protected:
  bool value_forceCompute;

  // filter parameters
  std::string m_operation;

  // filter output
  void *m_output_mesh_void;
};

} // namespace FEVV

