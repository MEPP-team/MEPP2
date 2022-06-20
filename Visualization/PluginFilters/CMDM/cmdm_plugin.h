// Copyright (c) 2012-2020 University of Lyon and CNRS (France).
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
#include "Dialogs/cmdm_dialog.h"

#ifndef Q_MOC_RUN // MT : very important to avoid the error : ' Parse error at
                  // "BOOST_JOIN" ' -> (qt4 pb with boost)
#include "Visualization/PluginFilters/BasePluginQt.h"
#include "Visualization/SimpleViewer.h"

#include "Visualization/SimpleWindow.h"

// A) include the header of the filter corresponding to your operation
#include "FEVV/Filters/CGAL/Surface_mesh/CMDM/cmdm.h"
#include "FEVV/Filters/Generic/color_mesh.h"
#include "FEVV/Filters/Generic/minmax_map.h"

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

class CMDMPlugin : public QObject,
                   public Generic_PluginInterface,
                   public BasePluginQt
{
  Q_OBJECT
  Q_INTERFACES(FEVV::Generic_PluginInterface)
#if(FEVV_USE_QT5) // see at the end of .cpp for QT4
  Q_PLUGIN_METADATA(IID "CMDMPlugin")
#endif

  /*public:
      using BasePlugin::apply;*/
public:
  CMDMPlugin() = default;
  ~CMDMPlugin() = default;

public:
  typedef typename FEVV::Vertex_pmap< MeshSurface, double > vertex_cmdm_map;

  using VertexColorMap =
      typename FEVV::PMap_traits< FEVV::vertex_color_t,
                                  FEVV::MeshSurface >::pmap_type;
  void init() override
  {
    one_two = false;
    two_one = false;
    scales = 3;
    color = false;
    LUT_CourbureClust = FEVV::Filters::make_LUT(false);
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

  void process(MeshSurface &m_degraded,
               FEVV::PMapsContainer &pmaps_bag_degraded,
               MeshSurface &m_original,
               FEVV::PMapsContainer &pmaps_bag_original,
               double &CMDM,
               vertex_cmdm_map &cmdm_pmap,
               VertexColorMap &v_cm)
  {
    auto pm_degrad = get(boost::vertex_point, m_degraded);
    auto pm_original = get(boost::vertex_point, m_original);

    using FaceNormalMap = typename FEVV::PMap_traits< FEVV::face_normal_t,
                                                      MeshSurface >::pmap_type;
    FaceNormalMap fnm_degrad;
    if(has_map(pmaps_bag_degraded, FEVV::face_normal))
    {
      std::cout << "use existing face-normal map for degraded mesh"
                << std::endl;
      fnm_degrad =
          get_property_map(FEVV::face_normal, m_degraded, pmaps_bag_degraded);
    }
    else
    {
      std::cout << "create face-normal map for degraded mesh" << std::endl;
      fnm_degrad = make_property_map(FEVV::face_normal, m_degraded);
      // store property map in property maps bag
      put_property_map(
          FEVV::face_normal, m_degraded, pmaps_bag_degraded, fnm_degrad);
      FEVV::Filters::calculate_face_normals(m_degraded, pm_degrad, fnm_degrad);
    }

    FaceNormalMap fnm_original;
    if(has_map(pmaps_bag_original, FEVV::face_normal))
    {
      std::cout << "use existing face-normal map for original mesh"
                << std::endl;
      fnm_original =
          get_property_map(FEVV::face_normal, m_original, pmaps_bag_original);
    }
    else
    {
      std::cout << "create face-normal map for original mesh" << std::endl;
      fnm_original = make_property_map(FEVV::face_normal, m_original);
      // store property map in property maps bag
      put_property_map(
          FEVV::face_normal, m_original, pmaps_bag_original, fnm_original);
      FEVV::Filters::calculate_face_normals(
          m_original, pm_original, fnm_original);
    }

    // Using a vertex property map to store the local distortions and visualize
    // them later (color map).
    // vertex_cmdm_map cmdm_pmap;
    if(has_map(pmaps_bag_degraded, std::string("v:cmdm")))
    {
      cmdm_pmap =
          boost::any_cast< vertex_cmdm_map >(pmaps_bag_degraded.at("v:cmdm"));
    }
    else
    {
      cmdm_pmap =
          FEVV::make_vertex_property_map< MeshSurface, double >(m_degraded);
      pmaps_bag_degraded["v:cmdm"] = cmdm_pmap;
    }

    // B) call the filter corresponding to your operation
    FEVV::Filters::process_CMDM_multires(m_degraded,
                                         pm_degrad,
                                         fnm_degrad,
                                         pmaps_bag_degraded,
                                         m_original,
                                         pm_original,
                                         fnm_original,
                                         pmaps_bag_original,
                                         scales,
                                         cmdm_pmap,
                                         CMDM);
    // To Visualise CMDM as color map later
    v_cm = get_property_map(FEVV::vertex_color, m_degraded, pmaps_bag_degraded);
  }

  template< typename HalfedgeGraph >
  void process(HalfedgeGraph &/*m_degraded*/,
               FEVV::PMapsContainer &/*pmaps_bag_degraded*/,
               HalfedgeGraph &/*m_original*/,
               FEVV::PMapsContainer &/*pmaps_bag_original*/,
               double &/*CMDM*/,
               vertex_cmdm_map &/*cmdm_pmap*/,
               VertexColorMap &/*v_cm*/)
  {
    QMessageBox::information(
        0, "", QObject::tr("CMDM only works with surface_mesh"));
    // return EXIT_FAILURE;
  }


  template< typename HalfedgeGraph >
  void applyHG(BaseAdapterVisu *_adapter,
               HalfedgeGraph * /*_mesh*/,
               FEVV::PMapsContainer * /*pmaps_bag*/)
  {
    // retrieve the two input meshes in current viewer window
    SimpleViewer *viewer =
        dynamic_cast< SimpleViewer * >(_adapter->getViewer());

    if(!viewer)
    {
      // something goes wrong
      QMessageBox::information(
          0, "", QObject::tr("An error occur during CMDM filter processing."));
      return;
    }

    MixedMeshesVector mixed_meshes = viewer->getMeshes();
    std::vector< FEVV::PMapsContainer * > properties_maps =
        viewer->get_properties_maps();

    if(mixed_meshes.size() == 2)
    {
      // check that the two input meshes have the same type
      if(mixed_meshes[0].second == mixed_meshes[1].second)
      {
        // get filter parameters from dialog window
        CMDMDialog dialog;
        if(dialog.exec() == QDialog::Accepted)
          dialog.getProcess(one_two, two_one, scales, color);
        else
          return; // abort applying filter

        // retieve meshes and property bags
        auto m1 = static_cast< HalfedgeGraph * >(mixed_meshes[0].first);
        auto pmaps_bag1 = properties_maps[0];
        auto m2 = static_cast< HalfedgeGraph * >(mixed_meshes[1].first);
        auto pmaps_bag2 = properties_maps[1];

        // ensure both meshes have diffuse color maps
        if(has_map(*pmaps_bag1, FEVV::vertex_color) &&
           has_map(*pmaps_bag2, FEVV::vertex_color))
        {
          // apply filter
          double CMDM = 0, CMDM_1_2 = 0, CMDM_2_1 = 0;
          vertex_cmdm_map cmdm_pmap_deg, cmdm_pmap_orig;
          VertexColorMap v_cm_deg, v_cm_orig;
          if(one_two)
          {
            process(*m2, *pmaps_bag2, *m1, *pmaps_bag1, CMDM_1_2, cmdm_pmap_deg,
                v_cm_deg);
            CMDM = CMDM_1_2;
          }
          if(two_one)
          {
            process(*m1, *pmaps_bag1, *m2, *pmaps_bag2, CMDM_2_1, cmdm_pmap_orig,
                v_cm_orig);
            CMDM = CMDM_2_1;
          }
          if(one_two && two_one)
          {
            CMDM = (CMDM_1_2 + CMDM_2_1) / 2.;
          }
          std::cout << "CMDM = " << CMDM << std::endl;

          // Visalize local distortions (LD)
          if(color)
          {
            double max_cmdm_1_2, min_cmdm_1_2, max_cmdm_2_1, min_cmdm_2_1;
            // 1- LD of the degraded mesh in relation to the original mesh(color
            // map of cmdm_1_2)
            FEVV::Filters::compute_min_max_vertices(*((MeshSurface *)(m2)),
                cmdm_pmap_deg,
                min_cmdm_1_2,
                max_cmdm_1_2);

            FEVV::Filters::color_vertices_from_map(*((MeshSurface *)(m2)),
                cmdm_pmap_deg,
                v_cm_deg,
                min_cmdm_1_2,
                max_cmdm_1_2,
                LUT_CourbureClust);

            // 2- LD of the original mesh in relation to the degraded mesh(color
            // map of cmdm_1_2)
            FEVV::Filters::compute_min_max_vertices(*((MeshSurface *)(m1)),
                cmdm_pmap_orig,
                min_cmdm_2_1,
                max_cmdm_2_1);

            FEVV::Filters::color_vertices_from_map(*((MeshSurface *)(m1)),
                cmdm_pmap_orig,
                v_cm_orig,
                min_cmdm_2_1,
                max_cmdm_2_1,
                LUT_CourbureClust);
          }

          // redraw meshes
          viewer->draw_or_redraw_mesh(m1, pmaps_bag1, true, false);
          viewer->draw_or_redraw_mesh(m2, pmaps_bag2, true, false);
        }
        else
        {
          QMessageBox::information(
              0,
              "",
              QObject::tr(
                  "CMDM needs 2 meshes with diffuse colors.\nUse MSDM2 Plugin"
                  "for meshes having geometry characteristics only."));
        }
      }
      else
      {
        QMessageBox::information(
            0,
            "",
            QObject::tr("CMDM filter can not be applied on meshes "
                        "with different datastructures."));
      }
    }
    else
    {
      QMessageBox::information(0,
                               "",
                               QObject::tr("CMDM needs 2 meshes opened "
                                           "with the same datastructure."));
    }

    reset();

#if(FEVV_USE_QT5)
    // empty
#else
    viewer->frame(); // necessary or not ?
#endif
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
    return QStringList() << "CMDMPlugin";
  }


  bool Generic_plugin(const QString &/*plugin*/) override
  {
    SimpleWindow *sw = static_cast< SimpleWindow * >(window);
    // dynamic_cast fails under OS X
    sw->onModificationParam("cmdm_qt_p", this);

    sw->onApplyButton();

    return true;
  }

signals:
  void resetSignal();

protected:
  bool one_two;
  bool two_one;
  int scales;
  bool color;
  FEVV::Filters::ColorMeshLUT LUT_CourbureClust;
};

} // namespace FEVV
