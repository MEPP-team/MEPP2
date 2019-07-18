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
#include "Visualization/PluginFilters/BasePluginQt.h"
#include "Visualization/SimpleViewer.h"

#include "Visualization/SimpleWindow.h"

#include "FEVV/Filters/CGAL/Boolean_Operations/boolean_operations.hpp"
#include "FEVV/Filters/Generic/homogeneous_transform.hpp"

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
                         public BasePluginQt
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
  void process(HalfedgeGraph *mesh_A,
               FEVV::PMapsContainer *pmaps_bag_A,
               Eigen::Matrix4d &matrix_A,
               HalfedgeGraph *mesh_B,
               FEVV::PMapsContainer *pmaps_bag_B,
               Eigen::Matrix4d &matrix_B)
  {
    std::cout << "Asking to apply BooleanOperations filter ! " << std::endl;

    // translate/rotate input meshes according to manipulators
    auto pm_A = get(boost::vertex_point, *mesh_A);
    FEVV::Filters::homogeneous_transform(*mesh_A, pm_A, matrix_A);
    auto pm_B = get(boost::vertex_point, *mesh_B);
    FEVV::Filters::homogeneous_transform(*mesh_B, pm_B, matrix_B);

    // normals of input meshes may be invalid due to mesh rotation
    remove_property_map(FEVV::vertex_normal, *pmaps_bag_A);
    remove_property_map(FEVV::face_normal, *pmaps_bag_A);
    remove_property_map(FEVV::vertex_normal, *pmaps_bag_B);
    remove_property_map(FEVV::face_normal, *pmaps_bag_B);

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

    // hide mesh A
    auto m_gpm_A =
        get_property_map(FEVV::mesh_guiproperties, *mesh_A, *pmaps_bag_A);
    auto gui_props_A = get(m_gpm_A, 0);
    //gui_props_A.is_visible = false; // we finally use TIME mode (see function 'activate_time_mode' below...)
    put(m_gpm_A, 0, gui_props_A);

    // hide mesh B
    auto m_gpm_B =
        get_property_map(FEVV::mesh_guiproperties, *mesh_B, *pmaps_bag_B);
    auto gui_props_B = get(m_gpm_B, 0);
    //gui_props_B.is_visible = false; // we finally use TIME mode (see function 'activate_time_mode' below...)
    put(m_gpm_B, 0, gui_props_B);

    // show output mesh
    FEVV::Types::GuiProperties gui_props_output;
    //gui_props_output.is_visible = true; // not necessary because true by default...

    // create a property map and a bag to store output mesh GUI properties
    auto m_gpm_output = make_property_map(FEVV::mesh_guiproperties, *output_mesh);
    put(m_gpm_output, 0, gui_props_output);
    output_pmaps_bag = new FEVV::PMapsContainer;
    put_property_map(FEVV::mesh_guiproperties,
                     *output_mesh,
                     *output_pmaps_bag,
                     m_gpm_output);
  }

  template< typename HalfedgeGraph >
  void applyHG(BaseAdapterVisu *_adapter,
               HalfedgeGraph *_mesh,
               FEVV::PMapsContainer *pmaps_bag)
  {
    // retrieve the two input meshes in current viewer window,
    // then apply the filter

    SimpleViewer *viewer =
      dynamic_cast< SimpleViewer * >(_adapter->getViewer());

    MixedMeshesVector mixed_meshes = viewer->getMeshes();
    std::vector< FEVV::PMapsContainer * > pmaps_bags =
        viewer->get_properties_maps();
    if(mixed_meshes.size() >= 2)
    {
      // check that the two input meshes have the same type
      if( mixed_meshes[0].second ==  mixed_meshes[1].second )
      {
        // get filter parameters from dialog window
        DialogBooleanOperations1 dial1;
        dial1.setParameters(m_operation);
        if(dial1.exec() == QDialog::Accepted)
          dial1.getParameters(m_operation);
        else
          return; // abort applying filter

        // retrieve the two input meshes
        auto mA = static_cast< HalfedgeGraph * >(mixed_meshes[0].first);
        auto pmaps_bagA = pmaps_bags[0];
        auto matrix44_A = viewer->getTransformMatrixEigen(0);

        auto mB = static_cast< HalfedgeGraph * >(mixed_meshes[1].first);
        auto pmaps_bagB = pmaps_bags[1];
        auto matrix44_B = viewer->getTransformMatrixEigen(1);

        // apply filter
        process(mA, pmaps_bagA, matrix44_A, mB, pmaps_bagB, matrix44_B);

        // reset transform matrix of A and B because transformation
        // is now applied to mesh coordinates
        viewer->resetTransformMatrix(0);
        viewer->resetTransformMatrix(1);
      }
      else
      {
        QMessageBox::information(0,
            "",
            QObject::tr("Boolean Operations filter "
                        "can not be applied on meshes "
                        "with different datastructures."));
      }
    }
    else
    {
      QMessageBox::information(0,
          "",
          QObject::tr("Boolean Operations filter needs "
                      "two meshes opened "
                      "with the same datastructure."));
    }

    // draw output mesh
    if(viewer)
    {
      // redraw input meshes
      // required because the user may have manually moved the meshes
      auto input_mesh_A = static_cast< HalfedgeGraph * >(mixed_meshes[0].first);
      auto input_mesh_B = static_cast< HalfedgeGraph * >(mixed_meshes[1].first);
      viewer->draw_or_redraw_mesh(input_mesh_A, pmaps_bags[0], true, false);
      viewer->draw_or_redraw_mesh(input_mesh_B, pmaps_bags[1], true, false);

      // draw output mesh
      auto output_mesh = static_cast< HalfedgeGraph * >( m_output_mesh_void);
      if(output_mesh)
      {
          // pmaps_bag is required for display
          viewer->draw_or_redraw_mesh(output_mesh,
                                      output_pmaps_bag,
                                      false,
                                      false,
                                      m_operation);
      }
    }

    //ELO comment next line to keep parameters between calls
    //reset();

    viewer->activate_time_mode();

    viewer->frame();
  }

#ifdef FEVV_USE_OPENMESH
#if 0 //TODO-elo-note  compiling error in Cartesian_converter.h with OM
  void apply(BaseAdapterVisu *_adapter,
             MeshOpenMesh *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    applyHG< MeshOpenMesh >(_adapter, _mesh, pmaps_bag);
  }
#endif
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
#if 0
	//TODO-elo-note: missing num_halfedges() used by CGAL::copy_face_graph()
  //               and compiling error in Cartesian_converter.h with AIF
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
    return QStringList() << "BooleanOperationsPlugin";
  }

  bool Generic_plugin(const QString &plugin) override
  {
    SimpleWindow *sw = static_cast< SimpleWindow * >(window);
      // dynamic_cast fails under OS X
    sw->onModificationParam("booleanoperations_qt_p", this);
    sw->onApplyButton();

    return true;
  }

signals:
  void resetSignal();

protected:
  bool value_forceCompute;

  // filter parameters
  std::string m_operation;

  // filter output
  void *m_output_mesh_void;
  FEVV::PMapsContainer *output_pmaps_bag;
};

} // namespace FEVV

