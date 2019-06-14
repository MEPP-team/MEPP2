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

#include "FEVV/DataStructures/DataStructures_cgal_point_set.h"


//
// SimpleViewer<>::internal_createMesh(...) specialization for FEVV::CGALPointSet
// with  HalfedgeGraph = FEVV::CGALPointSet and PointMap = FEVV::CGALPointSet
//
template< >
inline
void
FEVV::SimpleViewer::internal_createMesh< FEVV::CGALPointSet, FEVV::CGALPointSet >(
    osg::Geode *&geode,
    FEVV::CGALPointSet *_g,
    PMapsContainer *_pmaps,
    std::vector< osg::ref_ptr< osg::Geometry > >  &geometries,
    std::vector< osg::ref_ptr< osg::Geometry > >  &geometriesL,
    std::vector< osg::ref_ptr< osg::Geometry > >  &geometriesP,
    std::vector< osg::ref_ptr< osg::Geometry > >  &geometries_edges,
    std::vector< osg::ref_ptr< osg::Geometry > >  &geometries_vertices,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &vertexArrays,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &vertexArrays_edges,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &vertexArrays_vertices,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &normalsArrays,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &normalsArraysF,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &tangentsArrays,
    std::vector< osg::ref_ptr< osg::Vec4Array > > &colorsArrays,
    std::vector< osg::ref_ptr< osg::Vec4Array > > &colorsArrays_edges,
    std::vector< osg::ref_ptr< osg::Vec4Array > > &colorsArrays_vertices,
    std::vector< osg::ref_ptr< osg::Vec2Array > > &texcoordsArrays,
    FEVV::CGALPointSet *_pm,
    std::string _mesh_file)
{
  std::cout << "here FEVV::SimpleViewer< CGALPointSet >::internal_createMesh "
            << std::endl;

  using GraphTraits = boost::graph_traits< FEVV::CGALPointSet >;
  using GeometryTraits = FEVV::Geometry_traits< FEVV::CGALPointSet >;
  using vertex_descriptor = typename GraphTraits::vertex_descriptor;
  using vertex_iterator = typename GraphTraits::vertex_iterator;
  using Point = typename GeometryTraits::Point;

  //TODO-elo-rm  using GraphTraits = boost::graph_traits< FEVV::CGALPointSet >;
  //TODO-elo-rm  using GeometryTraits = FEVV::Geometry_traits< FEVV::CGALPointSet >;
  //TODO-elo-rm  using vertex_iterator = typename GraphTraits::vertex_iterator;
  //TODO-elo-rm  using vertex_descriptor = typename GraphTraits::vertex_descriptor;
  //TODO-elo-rm  using halfedge_point = typename GeometryTraits::Point;
  //TODO-elo-rm  using halfedge_vector = typename GeometryTraits::Vector;

  // property maps stuff
  using VertexNormalMap =
      typename FEVV::PMap_traits< FEVV::vertex_normal_t,
                                  FEVV::CGALPointSet >::pmap_type;
  using VertexColorMap = typename FEVV::PMap_traits< FEVV::vertex_color_t,
                                                     FEVV::CGALPointSet >::pmap_type;

  using MeshGuipropertiesMap =
      typename FEVV::PMap_traits< FEVV::mesh_guiproperties_t,
                                  FEVV::CGALPointSet >::pmap_type;

  VertexNormalMap v_nm;
  VertexColorMap v_cm;

  MeshGuipropertiesMap m_gpm;

  VertexNormalMap *_vt_nm = nullptr;
  VertexColorMap *_vt_cm = nullptr;

  MeshGuipropertiesMap *_m_gpm = nullptr;
  size_t _m_mm_size = 0;

  // retrieve or create mesh gui-properties property map
  if(has_map(*_pmaps, FEVV::mesh_guiproperties))
  {
    m_gpm = get_property_map(FEVV::mesh_guiproperties, *_g, *_pmaps);
    _m_gpm = &m_gpm;
  }
  else
  {
    m_gpm = make_property_map(FEVV::mesh_guiproperties, *_g);
    FEVV::Types::GuiProperties gui_props;
    put(m_gpm, 0, gui_props);
    put_property_map(FEVV::mesh_guiproperties, *_g, *_pmaps, m_gpm);
    _m_gpm = &m_gpm;
  }

  // --- vertex_normal
  if(has_map(*_pmaps, FEVV::vertex_normal))
  {
    if(!m_redraw)
      std::cout << "[SimpleViewer] using vertex normal" << std::endl;
    v_nm = get_property_map(FEVV::vertex_normal, *_g, *_pmaps);
    _vt_nm = &v_nm;
  }
#if 0 //TODO-elo-restore
  // compute vertex normals if not provided
  if(_vt_nm == nullptr || m_recomputeNT_if_redraw)
  {
    if(!m_redraw)
      std::cout << "[SimpleViewer] vertex normal missing, computing it"
                << std::endl;
    v_nm = make_property_map(FEVV::vertex_normal, *_g);
    _vt_nm = &v_nm;
    FEVV::Filters::calculate_vertex_normals(*_g, *_pm, f_nm, v_nm);
    put_property_map(
        FEVV::vertex_normal, *_g, *_pmaps, v_nm); // MT add (for later redraw)
  }
  // --- vertex_normal
#endif

  // --- vertex_color
  if(has_map(*_pmaps, FEVV::vertex_color))
  {
    if(!m_redraw)
      std::cout << "[SimpleViewer] using vertex color" << std::endl;

    v_cm = get_property_map(FEVV::vertex_color, *_g, *_pmaps);
    _vt_cm = &v_cm;

    if(!m_redraw)
    {
      m_UseVertexColor = true;
      m_UseFaceColor = false;
      m_UseTexture = false;
    }
  }
  // --- vertex_color

  if(m_redraw) // IMPORTANT -> ONLY IF REDRAW via GUI or CODE
  {
    if(m_UseVertexColor) // else automatic detection, as for a first DRAW
      if((!m_space_time) || (m_space_time && m_space_time_changeColorMode))
      {
        _vt_cm = nullptr;

        if(m_UseVertexColor)
        {
          if(has_map(*_pmaps, FEVV::vertex_color))
          {
            v_cm = get_property_map(FEVV::vertex_color, *_g, *_pmaps);
            _vt_cm = &v_cm;
          }
          /*else
          {
                  // -- TODO MT - TEMP - HERE WE FORCE a full RED vertex-color
          map for the REDRAW using Vector = typename GeometryTraits::Vector;

                  auto iterator_pair = vertices(*_g); // vertices() returns a
          vertex_iterator pair vertex_iterator vi = iterator_pair.first;
                  vertex_iterator vi_end = iterator_pair.second;
                  for (; vi != vi_end; ++vi)
                  {
                          v_cm[*vi] = Vector(1., 0., 0.); // RGB
                  }

                  _vt_cm = &v_cm;
                  // -- TODO MT - TEMP - HERE WE FORCE a full RED vertex-color
          map for the REDRAW
          }*/
          else
          {
            //_vt_nm = nullptr;
          }
        }
        else
        {
          //_vt_nm = nullptr;
          //_f_nm = nullptr;
        }
      }
  }

  QApplication::setOverrideCursor(Qt::ForbiddenCursor);
  SimpleWindow *sw = static_cast< SimpleWindow * >(
      getWindow()); // here static_cast instead of dynamic_cast only for OSX and
                    // because of plugins... don't understand why...
  sw->statusBar()->showMessage(
      QObject::tr("Create and populate OSG objects...") /*, 2000*/);

  GeometryTraits gt(*_g);

  if(geode == nullptr)
  {
    geode = new osg::Geode; // MT OSG_GEODE : bug from JL ? Ici geode devrait
                            // être osg::ref_ptr<osg::Geode> ???
  }
  else
  {
    geode->removeDrawables(0, geode->getNumDrawables());
  }

  // we must handle one geometry per texture/material
  if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
  {
    geometries.push_back(new osg::Geometry());
    geometries[0]->setSupportsDisplayList(true);
    geometries[0]->setUseDisplayList(true);
    geometries[0]->setUseVertexBufferObjects(
        false); /*geometries[0]->setUseVertexArrayObject(false);*/
  }
  if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
  {
    geometriesL.push_back(new osg::Geometry());
    geometriesL[0]->setSupportsDisplayList(true);
    geometriesL[0]->setUseDisplayList(true);
    geometriesL[0]->setUseVertexBufferObjects(
        false); /*geometriesL[0]->setUseVertexArrayObject(false);*/
  }
  if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
  {
    geometriesP.push_back(new osg::Geometry());
    geometriesP[0]->setSupportsDisplayList(true);
    geometriesP[0]->setUseDisplayList(true);
    geometriesP[0]->setUseVertexBufferObjects(
        false); /*geometriesP[0]->setUseVertexArrayObject(false);*/
  }
  if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
  {
    geometries_edges.push_back(new osg::Geometry());
    geometries_edges[0]->setSupportsDisplayList(true);
    geometries_edges[0]->setUseDisplayList(true);
    geometries_edges[0]->setUseVertexBufferObjects(
        false); /*geometries_edges[0]->setUseVertexArrayObject(false);*/
  }
  if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
  {
    geometries_vertices.push_back(new osg::Geometry());
    geometries_vertices[0]->setSupportsDisplayList(true);
    geometries_vertices[0]->setUseDisplayList(true);
    geometries_vertices[0]->setUseVertexBufferObjects(
        false); /*geometries_vertices[0]->setUseVertexArrayObject(false);*/
  }
  if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
  {
    vertexArrays.push_back(new osg::Vec3Array);
  }
  if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
  {
    vertexArrays_edges.push_back(new osg::Vec3Array);
  }
  if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
  {
    vertexArrays_vertices.push_back(new osg::Vec3Array);
  }
  if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
  {
    normalsArrays.push_back(new osg::Vec3Array);
  }
  if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
  {
    normalsArraysF.push_back(new osg::Vec3Array);
  }
  if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
  {
    tangentsArrays.push_back(new osg::Vec3Array);
  }
  if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
  {
    colorsArrays.push_back(new osg::Vec4Array);
  }
  if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
  {
    colorsArrays_edges.push_back(new osg::Vec4Array);
  }
  if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
  {
    colorsArrays_vertices.push_back(new osg::Vec4Array);
  }
  if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
  {
    texcoordsArrays.push_back(new osg::Vec2Array);
  }
  // std::cout << "---------> geometries.push_back" << std::endl;
  size_t mtl_id = 0;

  unsigned int sizeVertex = 0;
  unsigned int sizeSPoints = 0;

  Point p0;
  vertex_descriptor vd0;

  /// Adding vertices - superimpose and 'only_pts' mode only
  if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw)) // NEW
  {
    // std::cout << "---------> Adding vertices (for superimpose and 'only_pts'
    // mode)" << std::endl;

    VertexColorMap *SAVE_vt_cm = _vt_cm;

    if(has_map(*_pmaps, FEVV::vertex_color))
    {
      v_cm = get_property_map(FEVV::vertex_color, *_g, *_pmaps);
      _vt_cm = &v_cm;
    }

    // ---

    for(vertex_iterator v_it = vertices(*_g).first;
        v_it != vertices(*_g).second;
        ++v_it)
    {
      // retrieve vertex
      vd0 = *v_it;

      // retrieve vertex-point
      p0 = (*_pm)[vd0];
      vertexArrays_vertices[mtl_id]->push_back(
          Helpers::VectorConverter< FEVV::CGALPointSet >(p0));

      sizeSPoints++;

      // color
      if(_vt_cm)
        colorsArrays_vertices[mtl_id]->push_back(
            Helpers::VectorColorConverter< FEVV::CGALPointSet >(
                get(v_cm, *v_it))); // user/filter/plugin colors
      else
        colorsArrays_vertices[mtl_id]->push_back(
            Helpers::ColorConverter(Color::Green())); // default color
    }

    geometries_vertices[mtl_id]->addPrimitiveSet(new osg::DrawArrays(
        osg::PrimitiveSet::POINTS, 0, vertexArrays_vertices[mtl_id]->size()));

    _vt_cm = SAVE_vt_cm;
  }

  // NEW_HERE-01 (ADD)
  geometries_vertices[mtl_id]->setStateSet(NULL);
  {
    // set point size
    osg::ref_ptr< osg::Point > pt = new osg::Point();
    if(m_RenderSuperimposedVertices_Big)
      pt->setSize(5.0f);
    else if(m_RenderSuperimposedVertices)
      pt->setSize(3.0f);
    else
      pt->setSize(1.0f);
    geometries_vertices[mtl_id]->getOrCreateStateSet()->setAttribute(
        pt, osg::StateAttribute::ON);

    // light
    geometries_vertices[mtl_id]->getOrCreateStateSet()->setMode(
        GL_LIGHTING, osg::StateAttribute::OFF); // light always OFF for
                                                // superimpose vertices
  }

  sw->statusBar()->showMessage(QObject::tr("") /*, 2000*/);
  QApplication::restoreOverrideCursor();

  std::cout << "[SimpleViewer] I have drawn " << sizeVertex << " vertices." << std::endl;
  std::cout << "[SimpleViewer] I have also drawn " << sizeSPoints
            << " (superimpose) points." << std::endl;

  const auto loadingStartTime = std::chrono::system_clock::now();

  QApplication::setOverrideCursor(Qt::BusyCursor);
  sw->statusBar()->showMessage(QObject::tr("Render mesh...") /*, 2000*/);

  // NEW_HERE-02 (ADD)
  std::vector< osg::ref_ptr< osg::Vec3Array > > *_normalsArrays = nullptr;
  if(m_SmoothFlat_Shading)
    _normalsArrays = &normalsArrays;
  else
    _normalsArrays = &normalsArraysF;

  // NEW_HERE-03 (ADD)
  std::vector< osg::ref_ptr< osg::Geometry > > *_geometries = nullptr;
  _geometries = &geometriesP; // RenderMethod::RENDER_POINTS

  if(m_RenderMode == RenderMode::RENDER_SHADERS_DIRECT_LIGHTING ||
     m_RenderMode == RenderMode::RENDER_SHADERS_INDIRECT_LIGHTING)
  {
    using VertexTangentMap =
        typename FEVV::PMap_traits< FEVV::vertex_tangent_t,
                                    FEVV::CGALPointSet >::pmap_type;
    using FaceColorMap =
        typename FEVV::PMap_traits< FEVV::face_color_t,
                                    FEVV::CGALPointSet >::pmap_type;
    using VertexUVMap =
        typename FEVV::PMap_traits< FEVV::vertex_texcoord_t,
                                    FEVV::CGALPointSet >::pmap_type;
    using HalfedgeUVMap =
        typename FEVV::PMap_traits< FEVV::halfedge_texcoord_t,
                                    FEVV::CGALPointSet >::pmap_type;
    using MeshMaterialsMap =
        typename FEVV::PMap_traits< FEVV::mesh_materials_t,
                                    FEVV::CGALPointSet >::pmap_type;

    VertexTangentMap *v_tan_m = nullptr;
    FaceColorMap *_f_cm = nullptr;
    VertexUVMap *_vt_uv_m = nullptr;
    HalfedgeUVMap *_het_uv_m = nullptr;
    MeshMaterialsMap *_m_mm = nullptr;

    internal_loadShadedMesh(geode,
                            _g,
                            *_geometries,
                            geometries_edges,
                            geometries_vertices,
                            vertexArrays,
                            vertexArrays_edges,
                            vertexArrays_vertices,
                            *_normalsArrays,
                            tangentsArrays,
                            texcoordsArrays,
                            colorsArrays,
                            colorsArrays_edges,
                            colorsArrays_vertices,
                            _m_mm_size,
                            _vt_nm,
                            v_tan_m,
                            _vt_cm,
                            _f_cm,
                            _vt_uv_m,
                            _het_uv_m,
                            _m_mm);

    if(m_RenderMode == RenderMode::RENDER_SHADERS_INDIRECT_LIGHTING)
    {
      std::cout << "[SimpleViewer] Rendering with indirect lighting."
                << std::endl;
      // RM: sends uniform as true if we set the mode to indirect lighting
      geode->getOrCreateStateSet()->addUniform(
          new osg::Uniform("uniUseIndirectLighting", true));
    }
    else
    {
      std::cout << "[SimpleViewer] Rendering with direct lighting."
                << std::endl;
    }
  }
  else
  {
    using FaceColorMap =
        typename FEVV::PMap_traits< FEVV::face_color_t,
                                    FEVV::CGALPointSet >::pmap_type;
    using VertexUVMap =
        typename FEVV::PMap_traits< FEVV::vertex_texcoord_t,
                                    FEVV::CGALPointSet >::pmap_type;
    using HalfedgeUVMap =
        typename FEVV::PMap_traits< FEVV::halfedge_texcoord_t,
                                    FEVV::CGALPointSet >::pmap_type;
    using MeshMaterialsMap =
        typename FEVV::PMap_traits< FEVV::mesh_materials_t,
                                    FEVV::CGALPointSet >::pmap_type;

    FaceColorMap *_f_cm = nullptr;
    VertexUVMap *_vt_uv_m = nullptr;
    HalfedgeUVMap *_het_uv_m = nullptr;
    MeshMaterialsMap *_m_mm = nullptr;
    int _texture_type = NO_TEXCOORDS;

    internal_loadLegacyMesh(geode,
                            _g,
                            *_geometries,
                            geometries_edges,
                            geometries_vertices,
                            vertexArrays,
                            vertexArrays_edges,
                            vertexArrays_vertices,
                            *_normalsArrays,
                            texcoordsArrays,
                            colorsArrays,
                            colorsArrays_edges,
                            colorsArrays_vertices,
                            _m_mm_size,
                            _texture_type,
                            _vt_nm,
                            _vt_cm,
                            _f_cm,
                            _vt_uv_m,
                            _het_uv_m,
                            _m_mm);
  }

  sw->statusBar()->showMessage(QObject::tr("") /*, 2000*/);
  QApplication::restoreOverrideCursor();

  std::cout << "[SimpleViewer] Done 'loading' mesh (in graphic card) in "
            << std::chrono::duration_cast< std::chrono::duration< float > >(
                   std::chrono::system_clock::now() - loadingStartTime)
                   .count()
            << " seconds." << std::endl;

  gizmo->setNodeMask(m_ShowAxis ? 0xffffffff : 0x0);
  grid->setNodeMask(m_ShowGrid ? 0xffffffff : 0x0);

  BaseViewer *bv = dynamic_cast< BaseViewer * >(this);
  BaseAdapterVisuQt *bavQt =
      dynamic_cast< BaseAdapterVisuQt * >(bv->getAdapter());
  // geode->setName( std::string("Mesh ") +
  // std::to_string(Helpers::nbMeshDrawed++) + std::string(" [") +
  // bavQt->windowTitle().toStdString() + std::string("]") );

  std::string ds_name = FEVV::getDatastructureName(_g);
  geode->setName(_mesh_file + std::string(" [") + ds_name + " " +
                 bavQt->windowTitle().toStdString() + std::string("]"));
  geode->addDescription("MESH");
}

