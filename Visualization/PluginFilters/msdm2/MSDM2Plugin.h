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
#include "Dialogs/DialogMSDM21.h"

#ifndef Q_MOC_RUN // MT : very important to avoid the error : ' Parse error at
                  // "BOOST_JOIN" ' -> (qt4 pb with boost)
#include "Visualization/PluginFilters/BasePluginQt.h"
#include "Visualization/SimpleViewer.h"

#include "Visualization/SimpleWindow.h"

#include "FEVV/Filters/CGAL/Surface_mesh/msdm2.h" // A) include the header of the filter corresponding to your operation
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

class MSDM2Plugin : public QObject,
                    public Generic_PluginInterface,
                    public BasePluginQt
{
  Q_OBJECT
  Q_INTERFACES(FEVV::Generic_PluginInterface)
#if(FEVV_USE_QT5) // see at the end of .cpp for QT4
  Q_PLUGIN_METADATA(IID "MSDM2Plugin")
#endif

  /*public:
      using BasePlugin::apply;*/
public:
  MSDM2Plugin() = default;
  ~MSDM2Plugin() = default;

public:
  void init() override
  {
    one_two = false;
    two_one = false;
    scales = 3;
    color = true;
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
               double &MSDM2)
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


    typedef typename FEVV::Vertex_pmap< MeshSurface, double > vertex_msdm2_map;
    vertex_msdm2_map msdm2_pmap;

    if(has_map(pmaps_bag_degraded, std::string("v:msdm2")))
    {
      msdm2_pmap =
          boost::any_cast< vertex_msdm2_map >(pmaps_bag_degraded.at("v:msdm2"));
    }
    else
    {
      msdm2_pmap =
          FEVV::make_vertex_property_map< MeshSurface, double >(m_degraded);
      pmaps_bag_degraded["v:msdm2"] = msdm2_pmap;
    }

    FEVV::Filters::process_msdm2_multires(m_degraded,
                                          pm_degrad,
                                          fnm_degrad,
                                          pmaps_bag_degraded,
                                          m_original,
                                          pm_original,
                                          fnm_original,
                                          pmaps_bag_original,
                                          scales,
                                          msdm2_pmap,
                                          MSDM2);

    using VertexColorMap = typename FEVV::PMap_traits< FEVV::vertex_color_t,
                                                       MeshSurface >::pmap_type;
    VertexColorMap v_cm;
    if(has_map(pmaps_bag_degraded, FEVV::vertex_color))
    {
      std::cout << "use existing vertex-color map" << std::endl;
      v_cm =
          get_property_map(FEVV::vertex_color, m_degraded, pmaps_bag_degraded);
    }
    else
    {
      std::cout << "create vertex-color map" << std::endl;
      v_cm = make_property_map(FEVV::vertex_color, m_degraded);
      // store property map in property maps bag
      put_property_map(
          FEVV::vertex_color, m_degraded, pmaps_bag_degraded, v_cm);
    }

    double maxMSDM2, minMSDM2;
    if(color)
    {
      FEVV::Filters::compute_min_max_vertices(
          m_degraded, msdm2_pmap, minMSDM2, maxMSDM2);

      FEVV::Filters::color_vertices_from_map(
          m_degraded, msdm2_pmap, v_cm, minMSDM2, maxMSDM2, LUT_CourbureClust);
    }
  }

  template< typename HalfedgeGraph >
  void process(HalfedgeGraph &/*m_degraded*/,
               FEVV::PMapsContainer &/*pmaps_bag_degraded*/,
               HalfedgeGraph &/*m_original*/,
               FEVV::PMapsContainer &/*pmaps_bag_original*/,
               double &/*MSDM2*/)
  {
    QMessageBox::information(0,
        "",
        QObject::tr("MSMD2 filter only works with Surface_mesh "
                    "data structure"));
  }

  template< typename HalfedgeGraph >
  void applyHG(BaseAdapterVisu *_adapter,
               HalfedgeGraph * /*_mesh*/,
               FEVV::PMapsContainer * /*pmaps_bag*/)
  {
    // retrieve the two input meshes in current viewer window,
    // then apply the filter

    SimpleViewer *viewer =
        dynamic_cast< SimpleViewer * >(_adapter->getViewer());

    if(! viewer)
    {
      // something goes wrong
      QMessageBox::information(0,
          "",
          QObject::tr("An error occur during MSDM2 filter processing."));
      return;
    }

    MixedMeshesVector mixed_meshes = viewer->getMeshes();
    std::vector< FEVV::PMapsContainer * > properties_maps =
        viewer->get_properties_maps();

    if(mixed_meshes.size() == 2)
    {
      // check that the two input meshes have the same type
      if( mixed_meshes[0].second ==  mixed_meshes[1].second )
      {
        // get filter parameters from dialog window
        DialogMSDM21 dial1;
        if(dial1.exec() == QDialog::Accepted)
          dial1.getProcess(one_two, two_one, scales, color);
        else
          return; // abort applying filter

        auto m1 = static_cast< HalfedgeGraph * >(mixed_meshes[0].first);
        auto pm1 = properties_maps[0];
        auto m2 = static_cast< HalfedgeGraph * >(mixed_meshes[1].first);
        auto pm2 = properties_maps[1];

        double MSDM2 = 0, MSDM2_1_2 = 0, MSDM2_2_1 = 0;
        if(one_two)
        {
          process(*m2, *pm2, *m1, *pm1, MSDM2_1_2);
          MSDM2 = MSDM2_1_2;
        }
        if(two_one)
        {
          process(*m1, *pm1, *m2, *pm2, MSDM2_2_1);
          MSDM2 = MSDM2_2_1;
        }
        if(one_two && two_one)
        {
          MSDM2 = (MSDM2_1_2 + MSDM2_2_1) / 2.;
        }
        std::cout << "MSDM2 :" << MSDM2 << std::endl;

        // redraw meshes
        viewer->draw_or_redraw_mesh(m1, pm1, true, false);
        viewer->draw_or_redraw_mesh(m2, pm2, true, false);
      }
      else
      {
        QMessageBox::information(0,
            "",
            QObject::tr("MSDM2 filter can not be applied on meshes "
                        "with different data structures."));
      }
    }
    else
    {
      QMessageBox::information(0,
          "",
          QObject::tr("MSDM2 filter needs two meshes opened "
                      "with the Surface_mesh data structure."));
    }

    // comment next line to keep parameters between calls
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
    return QStringList() << "MSDM2Plugin";
  }

  bool Generic_plugin(const QString &/*plugin*/) override
  {
    SimpleWindow *sw = static_cast< SimpleWindow * >(window);
      // dynamic_cast fails under OS X
    sw->onModificationParam("msdm2_qt_p", this);
    sw->onApplyButton();

    return true;
  }

signals:
  void resetSignal();

protected:
  bool one_two;
  bool two_one;
  int  scales;
  bool color;
  FEVV::Filters::ColorMeshLUT LUT_CourbureClust;
};

} // namespace FEVV

