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
#include "Visualization/Helpers/OSGDebug.hpp"

//#include <osgFX/Outline>
//#include <osgFX/SpecularHighlights>
//#include <osg/LightModel>

// osgManipulator
#define MANIPULATOR
#include <osg/MatrixTransform>

#include <osgManipulator/TabBoxDragger>
//#include "Visualization/OSG/Manipulator/TabBoxDragger2.h"

#include <osgManipulator/TrackballDragger>
// osgManipulator

#include <osg/LineWidth>
#include <osg/Point>

#include <memory>

#include "Visualization/SimpleWindow.h"

#include "FEVV/Filters/Generic/calculate_face_normals.hpp"
#include "FEVV/Filters/Generic/Manifold/calculate_vertex_normals.hpp"
#include "FEVV/Filters/Generic/Manifold/calculate_vertices_tangent.hpp"
#include "FEVV/Filters/Generic/Manifold/calculate_halfedges_tangent.hpp"
#include "FEVV/Filters/Generic/translation.hpp"


inline
FEVV::SimpleViewer::SimpleViewer() : BaseViewerOSG()
{
#ifdef DEBUG_VISU2
  std::cout << "*** this=" << this << "    entering " << __func__ << std::endl;
#endif

#ifdef DEBUG_VISU2
  std::cout << "*** this=" << this << "    leaving " << __func__ << std::endl;
#endif
}


inline
FEVV::SimpleViewer::~SimpleViewer()
{
#ifdef DEBUG_VISU2
  std::cout << "*** this=" << this << "    entering " << __func__ << std::endl;
#endif

  // std::cout << "--> ~SimpleViewer" << std::endl;

  // BaseViewer *bv = dynamic_cast< BaseViewer* >(this);
  // BaseAdapterVisuQt *bavQt = dynamic_cast< BaseAdapterVisuQt*
  // >(bv->getAdapter());

  std::vector< Adapter * > *adapters = myWindow->getAdapters();
  for(unsigned int ii = 0; ii < (*adapters).size(); ++ii)
  {
    if((*adapters)[ii] == myAdapter)
    {
      (*adapters).erase((*adapters).begin() + ii);
      // std::cout << "--> myAdapter erase in ~SimpleViewer - TAB SIZE: " <<
      // (*adapters).size() << std::endl;

      SimpleWindow *sw = static_cast< SimpleWindow * >(getWindow());
      sw->update();
    }
  }

  /*for (unsigned int ii = 0; ii < v_geodes.size(); ++ii)
          delete v_geodes[ii]; // -> impossible, car "smart pointer"*/

  /*for (unsigned int ii = 0; ii < v_draggers1.size(); ++ii)
          delete v_draggers1[ii];	// -> impossible, car "smart pointer"
  for (unsigned int ii = 0; ii < v_draggers2.size(); ++ii)
          delete v_draggers2[ii];	// -> impossible, car "smart pointer"*/

  for(unsigned int ii = 0; ii < v_mixed_meshes.size(); ++ii)
  {
    {
      auto mesh_type_pair = v_mixed_meshes[ii];

#ifdef FEVV_USE_CGAL
      if(mesh_type_pair.second == "POLYHEDRON")
      {
        auto mesh_ptr = static_cast< FEVV::MeshPolyhedron* >(mesh_type_pair.first);
        std::cout << "[SimpleViewer] deleting mesh " << mesh_ptr  << " with datastructure POLYHEDRON" << std::endl;
        delete mesh_ptr;
      }
      if(mesh_type_pair.second == "SURFACEMESH")
      {
        auto mesh_ptr = static_cast< FEVV::MeshSurface* >(mesh_type_pair.first);
        std::cout << "[SimpleViewer] deleting mesh " << mesh_ptr  << " with datastructure SURFACEMESH" << std::endl;
        delete mesh_ptr;
      }
      if(mesh_type_pair.second == "LCC")
      {
        auto mesh_ptr = static_cast< FEVV::MeshLCC* >(mesh_type_pair.first);
        std::cout << "[SimpleViewer] deleting mesh " << mesh_ptr  << " with datastructure LCC" << std::endl;
        delete mesh_ptr;
      }
      if(mesh_type_pair.second == "CGALPOINTSET")
      {
        auto mesh_ptr = static_cast< FEVV::CGALPointSet* >(mesh_type_pair.first);
        std::cout << "[SimpleViewer] deleting mesh " << mesh_ptr
                  << " with datastructure CGALPOINTSET" << std::endl;
        delete mesh_ptr;
      }
#endif //FEVV_USE_CGAL

#ifdef FEVV_USE_OPENMESH
      if(mesh_type_pair.second == "OPENMESH")
      {
        auto mesh_ptr = static_cast< FEVV::MeshOpenMesh* >(mesh_type_pair.first);
        std::cout << "[SimpleViewer] deleting mesh " << mesh_ptr  << " with datastructure OPENMESH" << std::endl;
        delete mesh_ptr;
      }
#endif //FEVV_USE_OPENMESH

#ifdef FEVV_USE_AIF
      if(mesh_type_pair.second == "AIF")
      {
        auto mesh_ptr = static_cast< FEVV::MeshAIF* >(mesh_type_pair.first);
        std::cout << "[SimpleViewer] deleting mesh " << mesh_ptr  << " with datastructure AIF" << std::endl;
        delete mesh_ptr;
      }
#endif //FEVV_USE_AIF

#ifdef FEVV_USE_PCL
      if(mesh_type_pair.second == "PCLPOINTCLOUD")
      {
        auto mesh_ptr = static_cast< FEVV::PCLPointCloud* >(mesh_type_pair.first);
        std::cout << "[SimpleViewer] deleting mesh " << mesh_ptr
                  << " with datastructure PCLPOINTCLOUD" << std::endl;
        delete mesh_ptr;
      }
#endif //FEVV_USE_PCL
    }
  }
  
  for(unsigned int ii = 0; ii < v_properties_maps.size(); ++ii)
  {
    delete v_properties_maps[ii];
    // std::cout << "--> delete v_properties_maps[ii] (" <<
    // bavQt->windowTitle().toStdString() << ") in ~SimpleViewer" << std::endl;
  }

  /*for(unsigned int ii = 0; ii < v_geometries.size(); ++ii)
  {
    for(unsigned int jj = 0; jj < v_geometries[ii].size(); ++jj)
    {
      delete v_geometries[ii][jj]; // -> impossible, car "smart pointer"
      // std::cout << "--> delete v_geometries[ii][jj] (" <<
      // bavQt->windowTitle().toStdString() << ") in ~SimpleViewer" <<
  std::endl;
    }
  }*/

  // NOTE : same for geometries_edges, geometries_vertices, geometries_normals,
  //        vertexArrays, vertexArrays_edges, vertexArrays_vertices, vertexArrays_normals,
  //        colorsArrays, colorsArrays_edges, colorsArrays_vertices, colorsArrays_normals,
  //        etc...

#ifdef DEBUG_VISU2
  std::cout << "*** this=" << this << "    leaving " << __func__ << std::endl;
#endif
}


inline
void
FEVV::SimpleViewer::init()
{
  if(!Assert::check(
         !bIsInit, "is already init. Leaving...", "SimpleViewer::init"))
  {
    return;
  }

  // disable the default setting of viewer.done() by pressing Escape
  setKeyEventSetsDone(0);

  changeBackgroundColor(Color::WetAsphalt());

  // LIGHT (test)
#if(0)
  // osg:Light nous permet de donner à notre lumière ses caractéristiques
  osg::Light *pLight = new osg::Light;
  pLight->setLightNum(1); // ici cette lumière sera donc GL_LIGHT1

  // On va ici créer la lumière numéro 1 (GL_LIGHT1). C'est une lumière blanche
  // (couleur diffuse) et sa composante ambiante est également blanche.
  pLight->setAmbient(osg::Vec4d(1.0, 1.0, 1.0, 0.0));
  pLight->setDiffuse(osg::Vec4(1.0f, 1.0f, 1.0f, 0.0f));
  pLight->setSpecular(osg::Vec4(1.0f, 1.0f, 1.0f, 1.0f));

  // Elle est située à la position 2, 0, 2.
  // Le 4e paramètre w sert en fait ici, comme en OpenGL, à dire si notre
  // lumière est ponctuelle ou directionnelle. Si w=0 notre lumière est
  // directionnelle, sinon elle est ponctuelle. On a donc bien ici une lumière
  // ponctuelle située en 2,0,2.
  pLight->setPosition(osg::Vec4(2.0f, 0.0f, 2.0f, 1.0f));

  osg::LightSource *pLightSource = new osg::LightSource;
  pLightSource->setLight(pLight);
  root_node->addChild(pLightSource);

  // IMPORTANT : étant donné que le viewer donne une valeur par défaut à
  // GL_LIGHT0, nous allons éteindre celle-ci et activer la lumière numéro 1.
  osg::StateSet *state = root_node->getOrCreateStateSet();
  state->setMode(GL_LIGHT0, osg::StateAttribute::OFF);
  state->setMode(GL_LIGHT1, osg::StateAttribute::ON);
#endif
  // LIGHT (test)

  gizmo = Debug::createGizmo();
  gizmo->setNodeMask(m_ShowAxis ? 0xffffffff : 0x0);
  addGroup(gizmo);

  grid = Debug::createUnitGrid();
  grid->setNodeMask(m_ShowGrid ? 0xffffffff : 0x0);
  addGroup(grid);

  // hud = Debug::createHud(hudText);
  // addGroup(hud);

  // time
  i_time = current_i_time = -1;
  // time

  std::cout << "SimpleViewer init" << std::endl;

  bIsInit = true;
}


inline
bool
FEVV::SimpleViewer::isInit() const
{
  return bIsInit;
}


inline
bool
FEVV::SimpleViewer::isValid() const
{
  return bIsInit;
}

/*template< typename HalfedgeGraph >
bool
FEVV::SimpleViewer<HalfedgeGraph>::isSelected() const
{
    for( bool b : v_meshIsSelected )
    {
        if( b )
        {
            return true;
        }
    }
    return false;
}*/


inline
bool
FEVV::SimpleViewer::changeBackgroundColor(
    const FEVV::Color &_color)
{
  // osgViewer::View* _osgView = getView(0); // for osgViewer::CompositeViewer
  osgViewer::View *_osgView =
      dynamic_cast< osgViewer::View * >(this); // for osgViewer::Viewer
  _osgView->getCamera()->setClearColor(Helpers::ColorConverter(_color));

  return true;
}


inline
bool
FEVV::SimpleViewer::saveScreenshot(const std::string &_name)
{
  std::unique_ptr< osgViewer::ScreenCaptureHandler > scrn(
      new osgViewer::ScreenCaptureHandler());
  std::unique_ptr< osgViewer::ScreenCaptureHandler::WriteToFile > captureOper(
      new osgViewer::ScreenCaptureHandler::WriteToFile(_name, "png"));

  scrn->setCaptureOperation(captureOper.get());
  scrn->captureNextFrame(*this);
  this->frame();

  return true;
}


inline
void
FEVV::SimpleViewer::addModel(Model *_geode)
{
  root_node->addChild(_geode);
  if(myWindow != nullptr)
  {
    myWindow->notify();
  }

  // osgUtil::Optimizer optimizer;
  // optimizer.optimize( root_node );
}


inline
void
FEVV::SimpleViewer::addGroup(Group *_group)
{
  root_node->addChild(_group);
  if(myWindow != nullptr)
  {
    myWindow->notify();
  }

  // osgUtil::Optimizer optimizer;
  // optimizer.optimize( root_node );
}


inline
typename FEVV::SimpleViewer::DataModelVector *
FEVV::SimpleViewer::getDataModel()
{
  visitor->reset();

  root_node->traverse(*visitor);

  return visitor->exportResults();
}


inline
std::vector< osg::Geode * >
FEVV::SimpleViewer::getSelectedGeodes()
{
  std::vector< osg::Geode * > result;
  unsigned int i_pos = 0;
  for(bool b : v_meshIsSelected)
  {
    if(b)
    {
      result.push_back(v_geodes[i_pos]);
    }
    ++i_pos;
  }

  return result;
}


inline
std::vector< osg::Geode * >
FEVV::SimpleViewer::getGeodes()
{
  return v_geodes;
}


inline
FEVV::MixedMeshesVector
FEVV::SimpleViewer::getSelectedMeshes()
{
  FEVV::MixedMeshesVector result;
  unsigned int i_pos = 0;
  for(bool b : v_meshIsSelected)
  {
    if(b)
    {
      result.push_back(v_mixed_meshes[i_pos]);
    }
    ++i_pos;
  }

  return result;
}


inline
FEVV::MixedMeshesVector
FEVV::SimpleViewer::getMeshes()
{
  return v_mixed_meshes;
}


inline
std::vector< std::string >
FEVV::SimpleViewer::getSelectedMeshesNames()
{
  std::vector< std::string > result;
  unsigned int i_pos = 0;
  for(bool b : v_meshIsSelected)
  {
    if(b)
    {
      result.push_back(v_meshes_names[i_pos]);
    }
    ++i_pos;
  }

  return result;
}


inline
std::vector< std::string >
FEVV::SimpleViewer::getMeshesNames()
{
  return v_meshes_names;
}


inline
std::vector< FEVV::PMapsContainer * >
FEVV::SimpleViewer::getSelected_properties_maps()
{
  std::vector< FEVV::PMapsContainer * > result;
  unsigned int i_pos = 0;
  for(bool b : v_meshIsSelected)
  {
    if(b)
    {
      result.push_back(v_properties_maps[i_pos]);
    }
    ++i_pos;
  }

  return result;
}


inline
std::vector< FEVV::PMapsContainer * >
FEVV::SimpleViewer::get_properties_maps()
{
  return v_properties_maps;
}


inline
std::vector< osg::Group * >
FEVV::SimpleViewer::getSelectedDraggers1()
{
  std::vector< osg::Group * > result;
  unsigned int i_pos = 0;
  for(bool b : v_meshIsSelected)
  {
    if(b)
    {
      result.push_back(v_draggers1[i_pos]);
    }
    ++i_pos;
  }

  return result;
}


inline
std::vector< osg::Group * >
FEVV::SimpleViewer::getDraggers1()
{
  return v_draggers1;
}


inline
std::vector< osg::Group * >
FEVV::SimpleViewer::getSelectedDraggers2()
{
  std::vector< osg::Group * > result;
  unsigned int i_pos = 0;
  for(bool b : v_meshIsSelected)
  {
    if(b)
    {
      result.push_back(v_draggers2[i_pos]);
    }
    ++i_pos;
  }

  return result;
}


inline
std::vector< osg::Group * >
FEVV::SimpleViewer::getDraggers2()
{
  return v_draggers2;
}

#if 0 //TODO-elo-rm-?-ask_MTO
template< typename HalfedgeGraph >
HalfedgeGraph *
FEVV::SimpleViewer< HalfedgeGraph >::getMesh(unsigned int _position)
{
  Assert::check(_position < v_meshes.size(),
                "The given position must be lower than v_meshes.size()",
                "SimpleViewer::getMesh(int)");

  return v_meshes[_position];
}
#endif


template< typename HalfedgeGraph, typename PointMap >
osg::Geode *
FEVV::SimpleViewer::internal_createMesh(
    HalfedgeGraph *_g,
    PMapsContainer *_pmaps,
    std::vector< osg::ref_ptr< osg::Geometry > > &geometries,
    std::vector< osg::ref_ptr< osg::Geometry > > &geometriesL,
    std::vector< osg::ref_ptr< osg::Geometry > > &geometriesP,
    std::vector< osg::ref_ptr< osg::Geometry > > &geometries_edges,
    std::vector< osg::ref_ptr< osg::Geometry > > &geometries_vertices,
    std::vector< osg::ref_ptr< osg::Geometry > > &geometries_normals,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &vertexArrays,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &vertexArrays_edges,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &vertexArrays_vertices,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &vertexArrays_normals,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &normalsArrays,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &normalsArraysF,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &normalsArrays_edges,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &normalsArrays_vertices,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &tangentsArrays,
    std::vector< osg::ref_ptr< osg::Vec4Array > > &colorsArrays,
    std::vector< osg::ref_ptr< osg::Vec4Array > > &colorsArrays_edges,
    std::vector< osg::ref_ptr< osg::Vec4Array > > &colorsArrays_vertices,
    std::vector< osg::ref_ptr< osg::Vec4Array > > &colorsArrays_normals,
    std::vector< osg::ref_ptr< osg::Vec2Array > > &texcoordsArrays,
    PointMap *_pm,
    std::string _mesh_file)
{
  osg::Geode *geode =
      new osg::Geode; // MT OSG_GEODE : bug from JL ? Ici devrait être
                      // osg::ref_ptr<osg::Geode> geode = new osg::Geode(); ???
  internal_createMesh(geode,
                      _g,
                      _pmaps,
                      geometries,
                      geometriesL,
                      geometriesP,
                      geometries_edges,
                      geometries_vertices,
                      geometries_normals,
                      vertexArrays,
                      vertexArrays_edges,
                      vertexArrays_vertices,
                      vertexArrays_normals,
                      normalsArrays,
                      normalsArraysF,
                      normalsArrays_edges,
                      normalsArrays_vertices,
                      tangentsArrays,
                      colorsArrays,
                      colorsArrays_edges,
                      colorsArrays_vertices,
                      colorsArrays_normals,
                      texcoordsArrays,
                      _pm,
                      _mesh_file);
  return geode;
}


template< typename HalfedgeGraph, typename PointMap >
void
FEVV::SimpleViewer::internal_createMesh(
    osg::Geode *&geode,
    HalfedgeGraph *_g,
    PMapsContainer *_pmaps,
    std::vector< osg::ref_ptr< osg::Geometry > > &geometries,
    std::vector< osg::ref_ptr< osg::Geometry > > &geometriesL,
    std::vector< osg::ref_ptr< osg::Geometry > > &geometriesP,
    std::vector< osg::ref_ptr< osg::Geometry > > &geometries_edges,
    std::vector< osg::ref_ptr< osg::Geometry > > &geometries_vertices,
    std::vector< osg::ref_ptr< osg::Geometry > > &geometries_normals,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &vertexArrays,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &vertexArrays_edges,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &vertexArrays_vertices,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &vertexArrays_normals,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &normalsArrays,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &normalsArraysF,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &normalsArrays_edges,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &normalsArrays_vertices,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &tangentsArrays,
    std::vector< osg::ref_ptr< osg::Vec4Array > > &colorsArrays,
    std::vector< osg::ref_ptr< osg::Vec4Array > > &colorsArrays_edges,
    std::vector< osg::ref_ptr< osg::Vec4Array > > &colorsArrays_vertices,
    std::vector< osg::ref_ptr< osg::Vec4Array > > &colorsArrays_normals,
    std::vector< osg::ref_ptr< osg::Vec2Array > > &texcoordsArrays,
    PointMap *_pm,
    std::string _mesh_file)
{
  using GraphTraits = boost::graph_traits< HalfedgeGraph >;
  using GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph >;
  using face_iterator = typename GraphTraits::face_iterator;
  //using face_descriptor = typename GraphTraits::face_descriptor;
  using edge_iterator = typename GraphTraits::edge_iterator;
  //using edge_descriptor = typename GraphTraits::edge_descriptor;
  using halfedge_descriptor = typename GraphTraits::halfedge_descriptor;
  using vertex_iterator = typename GraphTraits::vertex_iterator;
  using vertex_descriptor = typename GraphTraits::vertex_descriptor;
  using halfedge_point = typename GeometryTraits::Point;
  using halfedge_vector = typename GeometryTraits::Vector;

  // property maps stuff
  using VertexNormalMap =
      typename FEVV::PMap_traits< FEVV::vertex_normal_t,
                                  HalfedgeGraph >::pmap_type;
  using FaceNormalMap = typename FEVV::PMap_traits< FEVV::face_normal_t,
                                                    HalfedgeGraph >::pmap_type;

  using VertexColorMap = typename FEVV::PMap_traits< FEVV::vertex_color_t,
                                                     HalfedgeGraph >::pmap_type;
  using FaceColorMap = typename FEVV::PMap_traits< FEVV::face_color_t,
                                                   HalfedgeGraph >::pmap_type;

  using VertexUVMap = typename FEVV::PMap_traits< FEVV::vertex_texcoord_t,
                                                  HalfedgeGraph >::pmap_type;
  using HalfedgeUVMap = typename FEVV::PMap_traits< FEVV::halfedge_texcoord_t,
                                                    HalfedgeGraph >::pmap_type;

  // RM: tangents are used primarily for normal mapping
  using VertexTangentMap =
      typename FEVV::PMap_traits< FEVV::vertex_tangent_t,
                                  HalfedgeGraph >::pmap_type;

  using FaceMaterialMap =
      typename FEVV::PMap_traits< FEVV::face_material_t,
                                  HalfedgeGraph >::pmap_type;
  using MeshMaterialsMap =
      typename FEVV::PMap_traits< FEVV::mesh_materials_t,
                                  HalfedgeGraph >::pmap_type;
  using MeshGuipropertiesMap =
      typename FEVV::PMap_traits< FEVV::mesh_guiproperties_t,
                                  HalfedgeGraph >::pmap_type;

  VertexNormalMap v_nm;
  FaceNormalMap f_nm;
  VertexColorMap v_cm;
  FaceColorMap f_cm;
  VertexUVMap v_uvm;
  HalfedgeUVMap h_uvm;

  VertexTangentMap vt_m; // RM

  FaceMaterialMap f_mm;
  MeshMaterialsMap m_mm;
  MeshGuipropertiesMap m_gpm;

  VertexNormalMap *_vt_nm = nullptr;
  FaceNormalMap *_f_nm = nullptr;
  VertexColorMap *_vt_cm = nullptr;
  FaceColorMap *_f_cm = nullptr;
  VertexUVMap *_vt_uv_m = nullptr;
  HalfedgeUVMap *_het_uv_m = nullptr;

  VertexTangentMap *v_tan_m = nullptr; // RM

  FaceMaterialMap *_f_mm = nullptr;
  MeshMaterialsMap *_m_mm = nullptr;
  MeshGuipropertiesMap *_m_gpm = nullptr;
  size_t _m_mm_size = 0;

  // textures stuff
  int _texture_type = NO_TEXCOORDS;

  FEVV::Filters::translate(*_g, *_pm, m_step, 0., 0.); // TEMP for test

  // retrieve or create mesh gui-properties property map
  if(has_map(*_pmaps, FEVV::mesh_guiproperties))
  {
    m_gpm = get_property_map(FEVV::mesh_guiproperties, *_g, *_pmaps);
    _m_gpm = &m_gpm;
    // std::cout << "[SimpleViewer] **************get mesh_guiproperties
    // property_map" << std::endl;
  }
  else
  {
    m_gpm = make_property_map(FEVV::mesh_guiproperties, *_g);
    FEVV::Types::GuiProperties gui_props;
    put(m_gpm, 0, gui_props);
    put_property_map(FEVV::mesh_guiproperties, *_g, *_pmaps, m_gpm);
    _m_gpm = &m_gpm;
    // std::cout << "[SimpleViewer] **************make mesh_guiproperties
    // property_map" << std::endl;
  }

  // --- face_normal
  if(has_map(*_pmaps, FEVV::face_normal))
  {
    if(!m_redraw)
      std::cout << "[SimpleViewer] using face normal" << std::endl;
    f_nm = get_property_map(FEVV::face_normal, *_g, *_pmaps);
    _f_nm = &f_nm;
  }
  // compute face normals if not provided
  if(_f_nm == nullptr || m_recomputeNT_if_redraw)
  {
    if(!m_redraw)
      std::cout << "[SimpleViewer] face normal missing, computing it"
                << std::endl;
    f_nm = make_property_map(FEVV::face_normal, *_g);
    _f_nm = &f_nm;
    FEVV::Filters::calculate_face_normals(*_g, *_pm, f_nm);
    put_property_map(
        FEVV::face_normal, *_g, *_pmaps, f_nm); // MT add (for later redraw)
  }
  // --- face_normal

  // --- vertex_normal
  if(has_map(*_pmaps, FEVV::vertex_normal))
  {
    if(!m_redraw)
      std::cout << "[SimpleViewer] using vertex normal" << std::endl;
    v_nm = get_property_map(FEVV::vertex_normal, *_g, *_pmaps);
    _vt_nm = &v_nm;
  }
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

  // --- face_color
  if(has_map(*_pmaps, FEVV::face_color))
  {
    if(!m_redraw)
      std::cout << "[SimpleViewer] using face color" << std::endl;
    f_cm = get_property_map(FEVV::face_color, *_g, *_pmaps);
    _f_cm = &f_cm;

    if(!m_redraw)
    {
      m_UseVertexColor = false;
      m_UseFaceColor = true;
      m_UseTexture = false;
    }
  }
  // --- face_color

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

  // --- texture
  if(has_map(*_pmaps, FEVV::vertex_texcoord))
  {
    if(!m_redraw)
      std::cout << "[SimpleViewer] using vertex texture-coordinates"
                << std::endl;
    v_uvm = get_property_map(FEVV::vertex_texcoord, *_g, *_pmaps);
    _vt_uv_m = &v_uvm;
    _texture_type = VERTEX_TEXCOORDS2D;

    if(!m_redraw)
    {
      m_UseVertexColor = false;
      m_UseFaceColor = false;
      m_UseTexture = true;
    }

    if(has_map(*_pmaps, FEVV::vertex_tangent) && (!m_recomputeNT_if_redraw))
    {
      if(!m_redraw)
        std::cout << "[SimpleViewer] using vertex tangents" << std::endl;
      vt_m = get_property_map(FEVV::vertex_tangent, *_g, *_pmaps);
      v_tan_m = &vt_m;
    }
    else
    {
      // RM: compute vertex tangents
      //   Note: shouldn't be added if no normal map available, causing extra
      //   process
      if(!m_redraw)
        std::cout << "[SimpleViewer] vertex tangents missing, computing it"
                  << std::endl;
      vt_m = make_property_map(FEVV::vertex_tangent, *_g);
      v_tan_m = &vt_m;
      FEVV::Filters::calculate_vertices_tangent(*_g, *_pm, v_uvm, *v_tan_m);
      put_property_map(FEVV::vertex_tangent, *_g, *_pmaps, *v_tan_m);
    }
  }
  else if(has_map(*_pmaps, FEVV::halfedge_texcoord))
  {
    if(!m_redraw)
      std::cout << "[SimpleViewer] using halfedge texture-coordinates"
                << std::endl;
    h_uvm = get_property_map(FEVV::halfedge_texcoord, *_g, *_pmaps);
    _het_uv_m = &h_uvm;
    _texture_type = HALFEDGE_TEXCOORDS2D;

    if(!m_redraw)
    {
      m_UseVertexColor = false;
      m_UseFaceColor = false;
      m_UseTexture = true;
    }

    if(has_map(*_pmaps, FEVV::vertex_tangent) && (!m_recomputeNT_if_redraw))
    {
      if(!m_redraw)
        std::cout << "[SimpleViewer] using halfedge tangents" << std::endl;
      vt_m = get_property_map(FEVV::vertex_tangent, *_g, *_pmaps);
      v_tan_m = &vt_m;
    }
    else
    {
      // RM: compute halfedge tangents
      //   Note: shouldn't be added if no normal map available, causing extra
      //   process
      if(!m_redraw)
        std::cout << "[SimpleViewer] halfedge tangents missing, computing it"
                  << std::endl;
      vt_m = make_property_map(FEVV::vertex_tangent, *_g);
      v_tan_m = &vt_m;
      FEVV::Filters::calculate_halfedges_tangent(*_g, *_pm, h_uvm, *v_tan_m);
      put_property_map(FEVV::vertex_tangent, *_g, *_pmaps, *v_tan_m);
    }
  }

  if(has_map(*_pmaps, FEVV::face_material))
  {
    f_mm = get_property_map(FEVV::face_material, *_g, *_pmaps);
    _f_mm = &f_mm;
  }
  else
  {
    if(_texture_type == HALFEDGE_TEXCOORDS2D)
      _texture_type =
          NO_TEXCOORDS; // if face material map is missing, we can not display
                        // any texture in HALFEDGE_TEXCOORDS2D mode
  }

  if(has_map(*_pmaps, FEVV::mesh_materials))
  {
    m_mm = get_property_map(FEVV::mesh_materials, *_g, *_pmaps);
    _m_mm = &m_mm;

    _m_mm_size = std::distance(m_mm.storage_begin(), m_mm.storage_end());
    if(!m_redraw)
      std::cout << "[SimpleViewer] number of materials: " << _m_mm_size
                << std::endl;
  }
  else
  {
    _texture_type = NO_TEXCOORDS; // disable textures if no material
    if(!m_redraw)
      std::cout << "[SimpleViewer] no material found, disabling textures"
                << std::endl;
  }
  // --- texture

  // TEMP - not used
  /*GeometryTraits gt(*_g);

  typename GeometryTraits::Point minAABB, maxAABB;

  FEVV::Tools::compute_bounding_box(*_g, *_pm, minAABB, maxAABB, gt);*/
  // TEMP - not used

  if(m_redraw) // IMPORTANT -> ONLY IF REDRAW via GUI or CODE
  {
    if(m_UseVertexColor || m_UseFaceColor ||
       m_UseTexture) // else automatic detection, as for a first DRAW
      if((!m_space_time) || (m_space_time && m_space_time_changeColorMode))
      {
        //_vt_nm = nullptr;
        //_f_nm = nullptr;
        _vt_cm = nullptr;
        _f_cm = nullptr;
        _vt_uv_m = nullptr;
        _het_uv_m = nullptr;

        // textures stuff
        _texture_type = NO_TEXCOORDS;

        _f_mm = nullptr;
        _m_mm = nullptr;
        //_m_mm_size = 0; // NEW
        // textures stuff

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
        else if(m_UseFaceColor)
        {
          if(has_map(*_pmaps, FEVV::face_color))
          {
            f_cm = get_property_map(FEVV::face_color, *_g, *_pmaps);
            _f_cm = &f_cm;
          }
          else
          {
            //_vt_nm = nullptr;
          }
        }
        else if(m_UseTexture)
        {
          if(has_map(*_pmaps, FEVV::vertex_texcoord))
          {
            v_uvm = get_property_map(FEVV::vertex_texcoord, *_g, *_pmaps);
            _vt_uv_m = &v_uvm;
            _texture_type = VERTEX_TEXCOORDS2D;
          }
          else if(has_map(*_pmaps, FEVV::halfedge_texcoord))
          {
            h_uvm = get_property_map(FEVV::halfedge_texcoord, *_g, *_pmaps);
            _het_uv_m = &h_uvm;
            _texture_type = HALFEDGE_TEXCOORDS2D;
          }
          else
          {
            //_vt_nm = nullptr;
          }

          if(has_map(*_pmaps, FEVV::face_material))
          {
            f_mm = get_property_map(FEVV::face_material, *_g, *_pmaps);
            _f_mm = &f_mm;
          }
          else
          {
            if(_texture_type == HALFEDGE_TEXCOORDS2D)
              _texture_type = NO_TEXCOORDS; // if face material map is missing,
                                            // we can not display any texture in
                                            // HALFEDGE_TEXCOORDS2D mode
          }

          if(has_map(*_pmaps, FEVV::mesh_materials))
          {
            m_mm = get_property_map(FEVV::mesh_materials, *_g, *_pmaps);
            _m_mm = &m_mm;

            _m_mm_size =
                std::distance(m_mm.storage_begin(), m_mm.storage_end());
          }
          else
            _texture_type = NO_TEXCOORDS; // disable textures if no material
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
    geometries_normals.push_back(new osg::Geometry());
    geometries_normals[0]->setSupportsDisplayList(true);
    geometries_normals[0]->setUseDisplayList(true);
    geometries_normals[0]->setUseVertexBufferObjects(
        false); /*geometries_normals[0]->setUseVertexArrayObject(false);*/
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
    vertexArrays_normals.push_back(new osg::Vec3Array);
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
    normalsArrays_edges.push_back(new osg::Vec3Array);
  }
  if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
  {
    normalsArrays_vertices.push_back(new osg::Vec3Array);
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
    colorsArrays_normals.push_back(new osg::Vec4Array);
  }
  if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
  {
    texcoordsArrays.push_back(new osg::Vec2Array);
  }
  // std::cout << "---------> geometries.push_back" << std::endl;
  size_t mtl_id = 0;

  // add ONLY if multi-textures
  for(size_t mi = 1; mi < _m_mm_size; mi++)
  {
    if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
    {
      geometries.push_back(new osg::Geometry());
      geometries[mi]->setSupportsDisplayList(true);
      geometries[mi]->setUseDisplayList(true);
      geometries[mi]->setUseVertexBufferObjects(
          false); /*geometries[mi]->setUseVertexArrayObject(false);*/
    }
    if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
    {
      geometriesL.push_back(new osg::Geometry());
      geometriesL[mi]->setSupportsDisplayList(true);
      geometriesL[mi]->setUseDisplayList(true);
      geometriesL[mi]->setUseVertexBufferObjects(
          false); /*geometriesL[mi]->setUseVertexArrayObject(false);*/
    }
    if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
    {
      geometriesP.push_back(new osg::Geometry());
      geometriesP[mi]->setSupportsDisplayList(true);
      geometriesP[mi]->setUseDisplayList(true);
      geometriesP[mi]->setUseVertexBufferObjects(
          false); /*geometriesP[mi]->setUseVertexArrayObject(false);*/
    }
    if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
    {
      geometries_edges.push_back(new osg::Geometry());
      geometries_edges[mi]->setSupportsDisplayList(true);
      geometries_edges[mi]->setUseDisplayList(true);
      geometries_edges[mi]->setUseVertexBufferObjects(
          false); /*geometries_edges[mi]->setUseVertexArrayObject(false);*/
    }
    if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
    {
      geometries_vertices.push_back(new osg::Geometry());
      geometries_vertices[mi]->setSupportsDisplayList(true);
      geometries_vertices[mi]->setUseDisplayList(true);
      geometries_vertices[mi]->setUseVertexBufferObjects(
          false); /*geometries_vertices[mi]->setUseVertexArrayObject(false);*/
    }
    if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
    {
      geometries_normals.push_back(new osg::Geometry());
      geometries_normals[mi]->setSupportsDisplayList(true);
      geometries_normals[mi]->setUseDisplayList(true);
      geometries_normals[mi]->setUseVertexBufferObjects(
          false); /*geometries_normals[mi]->setUseVertexArrayObject(false);*/
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
      vertexArrays_normals.push_back(new osg::Vec3Array);
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
      normalsArrays_edges.push_back(new osg::Vec3Array);
    }
    if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
    {
      normalsArrays_vertices.push_back(new osg::Vec3Array);
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
      colorsArrays_normals.push_back(new osg::Vec4Array);
    }
    if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
    {
      texcoordsArrays.push_back(new osg::Vec2Array);
    }
    // std::cout << "---------> geometries.push_back" << std::endl;
  }

  unsigned int sizeVertex = 0;
  unsigned int sizeFace = 0;
  unsigned int sizeSPoints = 0;
  unsigned int sizeSLines = 0;

  face_iterator fb, fe;
  halfedge_point p0;
  vertex_descriptor vd0;
  halfedge_vector normal;

  bool texture_corner_mode_on =
      (_het_uv_m != nullptr && _texture_type == HALFEDGE_TEXCOORDS2D);
  bool texture_vertex_mode_on =
      (_vt_uv_m != nullptr && _texture_type == VERTEX_TEXCOORDS2D);

  // std::cout << "---------> texture_corner_mode_on: " <<
  // texture_corner_mode_on << std::endl; std::cout << "--------->
  // texture_vertex_mode_on: " << texture_vertex_mode_on << std::endl;

  /// Adding edges - superimpose only
  // if(m_RenderSuperimposedEdges)
  if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw)) // NEW
  {
    const float MAGNITUDE = 0.0005;
    // std::cout << "---------> Adding edges (for superimpose)" << std::endl;

    using EdgeColorMap = typename FEVV::PMap_traits< FEVV::edge_color_t,
                                                     HalfedgeGraph >::pmap_type;
    EdgeColorMap e_cm;
    EdgeColorMap *_e_cm = nullptr;

    if(has_map(*_pmaps, FEVV::edge_color))
    {
      e_cm = get_property_map(FEVV::edge_color, *_g, *_pmaps);
      _e_cm = &e_cm;
    }

    // ---

    vertex_descriptor vs, vt;
    halfedge_point ps, pt;

    edge_iterator eb, ee;
    for(boost::tie(eb, ee) = edges(*_g); eb != ee; ++eb)
    {
      // retrieve vertex (source)
      vs = source(*eb, *_g);

      // retrieve vertex-point (source)
      ps = (*_pm)[vs];
      vertexArrays_edges[mtl_id]->push_back( Helpers::VectorConverter< HalfedgeGraph >(ps) + Helpers::VectorConverter< HalfedgeGraph >(get(v_nm, vs)) * MAGNITUDE ); // ok because _vt_nm always true

      // retrieve vertex (target)
      vt = target(*eb, *_g);

      // retrieve vertex-point (target)
      pt = (*_pm)[vt];
      vertexArrays_edges[mtl_id]->push_back( Helpers::VectorConverter< HalfedgeGraph >(pt) + Helpers::VectorConverter< HalfedgeGraph >(get(v_nm, vt)) * MAGNITUDE ); // ok because _vt_nm always true

      sizeSLines++;

      // color
      if(_e_cm)
      {
        colorsArrays_edges[mtl_id]->push_back(
            Helpers::VectorColorConverter< HalfedgeGraph >(
                get(e_cm, *eb))); // user/filter/plugin colors
        colorsArrays_edges[mtl_id]->push_back(
            Helpers::VectorColorConverter< HalfedgeGraph >(
                get(e_cm, *eb))); // user/filter/plugin colors
      }
      else
      {
        colorsArrays_edges[mtl_id]->push_back(
            Helpers::ColorConverter(Color::Yellow())); // default color
        colorsArrays_edges[mtl_id]->push_back(
            Helpers::ColorConverter(Color::Yellow())); // default color
      }

      // normal (only for GEOMETRY shader)
      if(_vt_nm)
      {
        normalsArrays_edges[mtl_id]->push_back(
            Helpers::VectorConverter< HalfedgeGraph >(
                get(v_nm, vs)));

        normalsArrays_edges[mtl_id]->push_back(
            Helpers::VectorConverter< HalfedgeGraph >(
                get(v_nm, vt)));
      }
    }

    geometries_edges[mtl_id]->addPrimitiveSet(new osg::DrawArrays(
        osg::PrimitiveSet::LINES, 0, vertexArrays_edges[mtl_id]->size()));

    // NEW_HERE-01 (DEL)
#if 0
    // set line width
    osg::ref_ptr< osg::LineWidth > linewidth = new osg::LineWidth();
    linewidth->setWidth(1.5f);
    geometries_edges[mtl_id]
      ->getOrCreateStateSet()
      ->setAttribute(linewidth, osg::StateAttribute::ON); // setAttributeAndModes (other function)

    // light
    geometries_edges[mtl_id]->getOrCreateStateSet()->setMode(
      GL_LIGHTING,
      osg::StateAttribute::OFF); // light always OFF for superimpose edges
#endif
  }
  // NEW_HERE-01 (ADD)
  geometries_edges[mtl_id]->setStateSet(NULL);
  {
    // set line width
    osg::ref_ptr< osg::LineWidth > linewidth = new osg::LineWidth();
    linewidth->setWidth(1.5f);
    geometries_edges[mtl_id]
        ->getOrCreateStateSet()
        ->setAttribute(linewidth, osg::StateAttribute::ON); // setAttributeAndModes (other function)

    // light
    geometries_edges[mtl_id]->getOrCreateStateSet()->setMode(
        GL_LIGHTING,
        osg::StateAttribute::OFF); // light always OFF for superimpose edges
  }

  /// Adding vertices - superimpose and 'only_pts' mode only
  size_t nb_faces = size_of_faces(*_g);
  // if(m_RenderSuperimposedVertices || m_RenderSuperimposedVertices_Big ||
  // (nb_faces==0)) // last test for 'only_pts' mode
  if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw)) // NEW
  {
    const float MAGNITUDE_N = 0.08;
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
          Helpers::VectorConverter< HalfedgeGraph >(p0));

      // [ normals
      vertexArrays_normals[mtl_id]->push_back( Helpers::VectorConverter< HalfedgeGraph >(p0) );
      vertexArrays_normals[mtl_id]->push_back( Helpers::VectorConverter< HalfedgeGraph >(p0) + Helpers::VectorConverter< HalfedgeGraph >(get(v_nm, *v_it)) * MAGNITUDE_N ); // ok because _vt_nm always true
      // ] normals

      sizeSPoints++;

      // color
      if(_vt_cm)
        colorsArrays_vertices[mtl_id]->push_back(
            Helpers::VectorColorConverter< HalfedgeGraph >(
                get(v_cm, *v_it))); // user/filter/plugin colors
      else
        colorsArrays_vertices[mtl_id]->push_back(
            Helpers::ColorConverter(Color::Green())); // default color

      // normal (only for GEOMETRY shader)
      if(_vt_nm)
        normalsArrays_vertices[mtl_id]->push_back(
            Helpers::VectorConverter< HalfedgeGraph >(
                get(v_nm, *v_it)));

      // [ normals
      colorsArrays_normals[mtl_id]->push_back(Helpers::ColorConverter(Color::Red()));
      colorsArrays_normals[mtl_id]->push_back(Helpers::ColorConverter(Color::Red()));
      // ] normals
    }

    geometries_vertices[mtl_id]->addPrimitiveSet(new osg::DrawArrays(
        osg::PrimitiveSet::POINTS, 0, vertexArrays_vertices[mtl_id]->size()));

    // [ normals
    geometries_normals[mtl_id]->addPrimitiveSet(new osg::DrawArrays(
        osg::PrimitiveSet::LINES, 0, vertexArrays_normals[mtl_id]->size()));
    // ] normals

    // NEW_HERE-01 (DEL)
#if 0
    if(m_RenderSuperimposedVertices || m_RenderSuperimposedVertices_Big)
    {
      // set point size
      osg::ref_ptr< osg::Point > pt = new osg::Point();
      if(m_RenderSuperimposedVertices_Big)
        pt->setSize(5.0f);
      else
        pt->setSize(3.0f);
      geometries_vertices[mtl_id]->getOrCreateStateSet()->setAttribute(
          pt, osg::StateAttribute::ON);
    }

    // light
    geometries_vertices[mtl_id]->getOrCreateStateSet()->setMode(
      GL_LIGHTING, osg::StateAttribute::OFF); // light always OFF for
                                              // superimpose vertices
#endif

    // ---

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
  // [ normals
  geometries_normals[mtl_id]->setStateSet(NULL);
  {
    // set line width
    osg::ref_ptr< osg::LineWidth > linewidth = new osg::LineWidth();
    linewidth->setWidth(1.5f);
    geometries_normals[mtl_id]
        ->getOrCreateStateSet()
        ->setAttribute(linewidth, osg::StateAttribute::ON); // setAttributeAndModes (other function)

    // light
    geometries_normals[mtl_id]->getOrCreateStateSet()->setMode(
        GL_LIGHTING,
        osg::StateAttribute::OFF); // light always OFF for normals
  }
  // ] normals

  /// Adding faces
  if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw)) // NEW
  {
    // std::cout << "---------> Adding faces" << std::endl;

    for(boost::tie(fb, fe) = faces(*_g); fb != fe; ++fb)
    {
      // retrieve the material/texture ID of the current face (ONLY in
      // HALFEDGE_TEXCOORDS2D mode)
      if(texture_corner_mode_on || texture_vertex_mode_on)
      {
        mtl_id = get(*_f_mm, *fb);
        // std::cout << "---------> mtl_id for current face: " << mtl_id <<
        // std::endl;
      }

      halfedge_descriptor edg = halfedge(*fb, *_g);
      halfedge_descriptor edg_begin = edg;

      std::vector< halfedge_point > p;
      std::vector< vertex_descriptor > vd;

      // loop over halfedges of the face
      do
      {
        vd.push_back(target(edg, *_g));
        p.push_back((*_pm)[vd.back()]);

        vertexArrays[mtl_id]->push_back(
            Helpers::VectorConverter< HalfedgeGraph >(p.back()));

        if(_vt_cm != nullptr)
        {
          colorsArrays[mtl_id]->push_back(
              Helpers::VectorColorConverter< HalfedgeGraph >(
                  (*_vt_cm)[vd.back()]));
        }

        if(_het_uv_m != nullptr && _texture_type == HALFEDGE_TEXCOORDS2D)
        {
          texcoordsArrays[mtl_id]->push_back(
              osg::Vec2((*_het_uv_m)[edg][0], (*_het_uv_m)[edg][1]));
        }
        else if(_vt_uv_m != nullptr && _texture_type == VERTEX_TEXCOORDS2D)
        {
          texcoordsArrays[mtl_id]->push_back(
              osg::Vec2((*_vt_uv_m)[vd.back()][0], (*_vt_uv_m)[vd.back()][1]));
        }

        if(_vt_nm != nullptr) // normal per vertex -> already calculated
        {
          // if(m_SmoothFlat_Shading) // NEW_HERE-02 (DEL)
          {
            normal = (*_vt_nm)[vd.back()];
            normalsArrays[mtl_id]->push_back(
                Helpers::VectorConverter< HalfedgeGraph >(
                    normal)); // smooth mode
          }
          // else // NEW_HERE-02 (DEL)
          {
            normal = (*_f_nm)[*fb];
            normalsArraysF[mtl_id]->push_back(
                Helpers::VectorConverter< HalfedgeGraph >(
                    normal)); // flat mode (caution, here, instead we take n x
                              // the normal per face)
          }

          if(v_tan_m != nullptr)
          {
            // RM: push tangents
            //   Note: shouldn't be added if no normal map available, causing
            //   extra process
            tangentsArrays[mtl_id]->push_back(
                osg::Vec3((*v_tan_m)[vd.back()][0],
                          (*v_tan_m)[vd.back()][1],
                          (*v_tan_m)[vd.back()][2]));
          }
        }

        ++sizeVertex;

        edg = next(edg, *_g);
      } while(edg != edg_begin);

      // create and populate OSG face object
      uint num_vertices_in_face = static_cast< uint >(p.size());

      auto drawing_method =
          static_cast< osg::PrimitiveSet::Mode >(m_RenderMethod);
      if(m_RenderMethod ==
         RenderMethod::RENDER_FILL) // RENDER_FILL is set to POLYGON primitive
      {
        if(num_vertices_in_face == 3)
          drawing_method = osg::PrimitiveSet::TRIANGLES;
        else if(num_vertices_in_face == 4)
        {
          drawing_method = osg::PrimitiveSet::QUADS;
          // std::cout << "---> QUAD-QUAD-QUAD-QUAD-QUAD" << std::endl;
        }
        else
        {
          // std::cout << "---> POLY-POLY-POLY-POLY-POLY" << std::endl;
        }
      }

      geometries[mtl_id]->addPrimitiveSet(new osg::DrawArrays(
          drawing_method,
          vertexArrays[mtl_id]->size() - num_vertices_in_face,
          num_vertices_in_face));
      // NEW_HERE-03 (ADD)
      geometriesL[mtl_id]->addPrimitiveSet(new osg::DrawArrays(
          RenderMethod::RENDER_LINES,
          vertexArrays[mtl_id]->size() - num_vertices_in_face,
          num_vertices_in_face));
      geometriesP[mtl_id]->addPrimitiveSet(new osg::DrawArrays(
          RenderMethod::RENDER_POINTS,
          vertexArrays[mtl_id]->size() - num_vertices_in_face,
          num_vertices_in_face));

      // populate face normal array
      if(_vt_nm !=
         nullptr) // normal per vertex (see above) -> already calculated
      {
        // do nothing
        // std::cout << "---> normal per vertex - normal per vertex - normal per
        // vertex" << std::endl;
      }
      else if(_f_nm != nullptr) // normal per face -> already calculated
      {
        std::cout << "---> normal per face - normal per face - normal per face --> MUST NEVER HAPPENS !!!"
                  << std::endl; // TEMP

        normal = (*_f_nm)[*fb];
        normalsArrays[mtl_id]->push_back(
            Helpers::VectorConverter< HalfedgeGraph >(normal));
      }
      else // if none - BUT NEVER HAPPENS - re-compute normal for current face
      {
        std::cout << "---> NEVER HAPPENS - NEVER HAPPENS - NEVER HAPPENS"
                  << std::endl; // TEMP

#if 0
      if(num_vertices_in_face == 3) // TRIANGLE
      {
        normal = gt.unit_normal(p[1], p[2], p[0]);
        normalsArrays[mtl_id]->push_back(
            Helpers::VectorConverter< HalfedgeGraph >(normal));
      }
      else if(num_vertices_in_face == 4) // QUAD
      {
        normal = gt.NULL_VECTOR;
        normal = gt.normal(p[1], p[2], p[0]) + gt.normal(p[2], p[3], p[1]);
        normal = normal / gt.length(normal);
        normalsArrays[mtl_id]->push_back(
            Helpers::VectorConverter< HalfedgeGraph >(normal));
      }
      else // POLYGON
      {
        // TODO
        std::cerr << "-----> [SimpleViewer] Only re-computed normals for "
                     "triangles/quads are implemented -> re-computed normals "
                     "for other polygons are not yet implemented."
                  << std::endl;
      }
#endif
      }

      // populate face color array
      if(_f_cm != nullptr)
      {
        colorsArrays[mtl_id]->push_back(
            Helpers::VectorColorConverter< HalfedgeGraph >((*_f_cm)[*fb]));
      }

      ++sizeFace;
    }

    std::cout << "[SimpleViewer] I have drawn " << sizeFace << " faces (with "
            << sizeVertex << " vertices)." << std::endl;
    std::cout << "[SimpleViewer] I have also drawn " << sizeSPoints
            << " (superimpose) points and " << sizeSLines
            << " superimpose lines." << std::endl;
  }

  sw->statusBar()->showMessage(QObject::tr("") /*, 2000*/);
  QApplication::restoreOverrideCursor();

  // auto gui_props = get((*_m_gpm), 0);
  // geode->setNodeMask(gui_props.is_visible ? 0xffffffff : 0x0); // 19/03/19

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
  if(m_RenderMethod == RenderMethod::RENDER_FILL)
    _geometries = &geometries;
  else if(m_RenderMethod == RenderMethod::RENDER_LINES)
    _geometries = &geometriesL;
  else
    _geometries = &geometriesP; // RenderMethod::RENDER_POINTS

  if(m_RenderMode == RenderMode::RENDER_SHADERS_DIRECT_LIGHTING ||
     m_RenderMode == RenderMode::RENDER_SHADERS_INDIRECT_LIGHTING)
  {
    internal_loadShadedMesh(geode,
                            _g,
                            *_geometries,
                            geometries_edges,
                            geometries_vertices,
                            geometries_normals,
                            vertexArrays,
                            vertexArrays_edges,
                            vertexArrays_vertices,
                            vertexArrays_normals,
                            *_normalsArrays,
                            normalsArrays_edges,
                            normalsArrays_vertices,
                            tangentsArrays,
                            texcoordsArrays,
                            colorsArrays,
                            colorsArrays_edges,
                            colorsArrays_vertices,
                            colorsArrays_normals,
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
    internal_loadLegacyMesh(geode,
                            _g,
                            *_geometries,
                            geometries_edges,
                            geometries_vertices,
                            geometries_normals,
                            vertexArrays,
                            vertexArrays_edges,
                            vertexArrays_vertices,
                            vertexArrays_normals,
                            *_normalsArrays,
                            normalsArrays_edges,
                            normalsArrays_vertices,
                            texcoordsArrays,
                            colorsArrays,
                            colorsArrays_edges,
                            colorsArrays_vertices,
                            colorsArrays_normals,
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
                 bavQt->windowTitle().left(bavQt->windowTitle().indexOf('>') + 1).toStdString() + std::string("]"));
  geode->addDescription("MESH");
}


template< typename PointCloud, typename PointMap >
void
FEVV::SimpleViewer::internal_createMesh_pointcloud(
    osg::Geode *&geode,
    PointCloud *_g,
    PMapsContainer *_pmaps,
    std::vector< osg::ref_ptr< osg::Geometry > > &geometries,
    std::vector< osg::ref_ptr< osg::Geometry > > &geometriesL,
    std::vector< osg::ref_ptr< osg::Geometry > > &geometriesP,
    std::vector< osg::ref_ptr< osg::Geometry > > &geometries_edges,
    std::vector< osg::ref_ptr< osg::Geometry > > &geometries_vertices,
    std::vector< osg::ref_ptr< osg::Geometry > > &geometries_normals,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &vertexArrays,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &vertexArrays_edges,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &vertexArrays_vertices,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &vertexArrays_normals,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &normalsArrays,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &normalsArraysF,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &normalsArrays_edges,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &normalsArrays_vertices,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &tangentsArrays,
    std::vector< osg::ref_ptr< osg::Vec4Array > > &colorsArrays,
    std::vector< osg::ref_ptr< osg::Vec4Array > > &colorsArrays_edges,
    std::vector< osg::ref_ptr< osg::Vec4Array > > &colorsArrays_vertices,
    std::vector< osg::ref_ptr< osg::Vec4Array > > &colorsArrays_normals,
    std::vector< osg::ref_ptr< osg::Vec2Array > > &texcoordsArrays,
    PointMap *_pm,
    std::string _mesh_file)
{
  using GraphTraits = boost::graph_traits< PointCloud >;
  using GeometryTraits = FEVV::Geometry_traits< PointCloud >;
  using vertex_descriptor = typename GraphTraits::vertex_descriptor;
  using vertex_iterator = typename GraphTraits::vertex_iterator;
  using Point = typename GeometryTraits::Point;
  using Vector = typename GeometryTraits::Vector;

  // property maps stuff
  using VertexNormalMap =
      typename FEVV::PMap_traits< FEVV::vertex_normal_t,
                                  PointCloud >::pmap_type;
  using VertexColorMap = typename FEVV::PMap_traits< FEVV::vertex_color_t,
                                                     PointCloud >::pmap_type;

  using MeshGuipropertiesMap =
      typename FEVV::PMap_traits< FEVV::mesh_guiproperties_t,
                                  PointCloud >::pmap_type;

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
    geometries_normals.push_back(new osg::Geometry());
    geometries_normals[0]->setSupportsDisplayList(true);
    geometries_normals[0]->setUseDisplayList(true);
    geometries_normals[0]->setUseVertexBufferObjects(
        false); /*geometries_normals[0]->setUseVertexArrayObject(false);*/
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
    vertexArrays_normals.push_back(new osg::Vec3Array);
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
    normalsArrays_edges.push_back(new osg::Vec3Array);
  }
  if((!m_redraw) || (m_redraw && m_recreateOSGobj_if_redraw))
  {
    normalsArrays_vertices.push_back(new osg::Vec3Array);
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
    colorsArrays_normals.push_back(new osg::Vec4Array);
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
    const float MAGNITUDE_N = 0.08;
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
          Helpers::VectorConverter< PointCloud >(
              Vector(GeometryTraits::get_x(p0),
                     GeometryTraits::get_y(p0),
                     GeometryTraits::get_z(p0))));
      // [ normals
      vertexArrays_normals[mtl_id]->push_back( Helpers::VectorConverter< PointCloud >(p0) );
      vertexArrays_normals[mtl_id]->push_back( Helpers::VectorConverter< PointCloud >(p0) + Helpers::VectorConverter< PointCloud >(get(v_nm, *v_it)) * MAGNITUDE_N ); // ok because _vt_nm always true
      // ] normals

      sizeSPoints++;

      // color
      if(_vt_cm)
      {
        auto color_uint8 = get(v_cm, *v_it);
        osg::Vec4 colorf((float)(color_uint8[0]/255.0f),
                         (float)(color_uint8[1]/255.0f),
                         (float)(color_uint8[2]/255.0f),
                         1.0f);
        colorsArrays_vertices[mtl_id]->push_back(colorf);
      }
      else
      {
        colorsArrays_vertices[mtl_id]->push_back(
            Helpers::ColorConverter(Color::Green())); // default color
      }

      // normal (only for GEOMETRY shader)
      if(_vt_nm)
        normalsArrays_vertices[mtl_id]->push_back(
            Helpers::VectorConverter< PointCloud >(
                get(v_nm, *v_it)));

      // [ normals
      colorsArrays_normals[mtl_id]->push_back(Helpers::ColorConverter(Color::Red()));
      colorsArrays_normals[mtl_id]->push_back(Helpers::ColorConverter(Color::Red()));
      // ] normals
    }

    geometries_vertices[mtl_id]->addPrimitiveSet(new osg::DrawArrays(
        osg::PrimitiveSet::POINTS, 0, vertexArrays_vertices[mtl_id]->size()));

    // [ normals
    geometries_normals[mtl_id]->addPrimitiveSet(new osg::DrawArrays(
        osg::PrimitiveSet::LINES, 0, vertexArrays_normals[mtl_id]->size()));
    // ] normals

    _vt_cm = SAVE_vt_cm;

    std::cout << "[SimpleViewer] I have drawn " << sizeSPoints << " points." << std::endl;
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

  // [ normals
  geometries_normals[mtl_id]->setStateSet(NULL);
  {
    // set line width
    osg::ref_ptr< osg::LineWidth > linewidth = new osg::LineWidth();
    linewidth->setWidth(1.5f);
    geometries_normals[mtl_id]
        ->getOrCreateStateSet()
        ->setAttribute(linewidth, osg::StateAttribute::ON); // setAttributeAndModes (other function)

    // light
    geometries_normals[mtl_id]->getOrCreateStateSet()->setMode(
        GL_LIGHTING,
        osg::StateAttribute::OFF); // light always OFF for normals
  }
  // ] normals

  sw->statusBar()->showMessage(QObject::tr("") /*, 2000*/);
  QApplication::restoreOverrideCursor();

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
                                    PointCloud >::pmap_type;
    using FaceColorMap =
        typename FEVV::PMap_traits< FEVV::face_color_t,
                                    PointCloud >::pmap_type;
    using VertexUVMap =
        typename FEVV::PMap_traits< FEVV::vertex_texcoord_t,
                                    PointCloud >::pmap_type;
    using HalfedgeUVMap =
        typename FEVV::PMap_traits< FEVV::halfedge_texcoord_t,
                                    PointCloud >::pmap_type;
    using MeshMaterialsMap =
        typename FEVV::PMap_traits< FEVV::mesh_materials_t,
                                    PointCloud >::pmap_type;

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
                            geometries_normals,
                            vertexArrays,
                            vertexArrays_edges,
                            vertexArrays_vertices,
                            vertexArrays_normals,
                            *_normalsArrays,
                            normalsArrays_edges,
                            normalsArrays_vertices,
                            tangentsArrays,
                            texcoordsArrays,
                            colorsArrays,
                            colorsArrays_edges,
                            colorsArrays_vertices,
                            colorsArrays_normals,
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
                                    PointCloud >::pmap_type;
    using VertexUVMap =
        typename FEVV::PMap_traits< FEVV::vertex_texcoord_t,
                                    PointCloud >::pmap_type;
    using HalfedgeUVMap =
        typename FEVV::PMap_traits< FEVV::halfedge_texcoord_t,
                                    PointCloud >::pmap_type;
    using MeshMaterialsMap =
        typename FEVV::PMap_traits< FEVV::mesh_materials_t,
                                    PointCloud >::pmap_type;

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
                            geometries_normals,
                            vertexArrays,
                            vertexArrays_edges,
                            vertexArrays_vertices,
                            vertexArrays_normals,
                            *_normalsArrays,
                            normalsArrays_edges,
                            normalsArrays_vertices,
                            texcoordsArrays,
                            colorsArrays,
                            colorsArrays_edges,
                            colorsArrays_vertices,
                            colorsArrays_normals,
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
                 bavQt->windowTitle().left(bavQt->windowTitle().indexOf('>') + 1).toStdString() + std::string("]"));
  geode->addDescription("MESH");
}


// osgManipulator
class ScaleConstraint : public osgManipulator::Constraint
{
public:
  ScaleConstraint() {}

  virtual bool constrain(osgManipulator::Scale1DCommand &command) const
  {
    command.setScale(1.0f);
    // OSG_NOTICE << "ScaleConstraint Scale1DCommand" << command.getScale() <<
    // std::endl;
    return true;
  }
  virtual bool constrain(osgManipulator::Scale2DCommand &command) const
  {
    command.setScale(osg::Vec2d(1.0, 1.0));
    // OSG_NOTICE << "ScaleConstraint Scale2DCommand " << command.getScale() <<
    // std::endl;
    return true;
  }
};


inline
osgManipulator::Dragger *
createDragger(const std::string &name)
{
  osgManipulator::Dragger *dragger = 0;

  if(name == "TrackballDragger")
  {
    osgManipulator::TrackballDragger *d =
        new osgManipulator::TrackballDragger();
    d->setupDefaultGeometry();
#if OSG_MIN_VERSION_REQUIRED(3, 4, 0)
    d->setAxisLineWidth(2.0f); // commented by default
                               // d->setPickCylinderHeight(0.1f);
#endif
    dragger = d;
  }
  else // TabBoxDragger (by default)
  {
    osgManipulator::TabBoxDragger *d = new osgManipulator::TabBoxDragger();
    // osgManipulator::TabBoxDragger2* d = new osgManipulator::TabBoxDragger2();
    d->setupDefaultGeometry();
    d->setPlaneColor(osg::Vec4(1.0f, 1.0f, 0.0f, 1.0f)); // yellow
    d->addConstraint(new ScaleConstraint());
    dragger = d;
  }

  dragger->setName(name);

  return dragger;
}


inline
osg::Node *
FEVV::SimpleViewer::addDraggersToScene(
    osg::Node *scene,
    const std::string &nameDrag1,
    float fScaleDrag1,
    char keyDrag1,
    const std::string &nameDrag2,
    float fScaleDrag2,
    char keyDrag2)
{
  // scene->getOrCreateStateSet()->setMode(GL_NORMALIZE,
  // osg::StateAttribute::ON); // necessary ? // this ensures correct lighting
  // for scaled draggers


  // osg::Group "GroupRoot"
  //    |
  //    |__ osg::MatrixTransform "MatrixTransform"
  //    |   |
  //    |   |__ osg::Group "DraggerGrp2_rotate"
  //    |       |
  //    |       |__ osgManipulator::TrackballDragger "TrackballDragger"
  //    |
  //    |__ osg::Group "DraggerGrp1_translate"
  //        |
  //        |__ osgManipulator::TabBoxDragger "TabBoxDragger"


  osg::MatrixTransform *transform = new osg::MatrixTransform;
  transform->setName("MatrixTransform");
  transform->addChild(scene);

  osg::Group *root = new osg::Group;
  root->setName("GroupRoot");
  root->addChild(transform);

  float scale;

  osgManipulator::Dragger *dragger1 = createDragger(nameDrag1);
  // here, we create a group only for the DataVisitor
  osg::Group *draggergroup1 = new osg::Group;
  draggergroup1->setName("DraggerGrp1_translate");
  draggergroup1->addChild(dragger1);
  root->addChild(draggergroup1);
  scale = scene->getBound().radius() * fScaleDrag1;
  dragger1->setMatrix(osg::Matrix::scale(scale, scale, scale) *
                      osg::Matrix::translate(scene->getBound().center()));
  dragger1->addTransformUpdating(transform);
  dragger1->setHandleEvents(true);
  dragger1->setActivationKeyEvent(keyDrag1);
  // dragger1->setDraggerActive(false); // not working ?
  dragger1->setNodeMask(m_ShowTranslateDragger ? 0xffffffff : 0x0);
  v_draggers1.push_back(dragger1);

  osgManipulator::Dragger *dragger2 = createDragger(nameDrag2);
  // here, we create a group only for the DataVisitor
  osg::Group *draggergroup2 = new osg::Group;
  draggergroup2->setName("DraggerGrp2_rotate");
  draggergroup2->addChild(dragger2);
  transform->addChild(
      draggergroup2); // IMPORTANT for the SECOND dragger -> NOT addChild from
                      // 'root' here BUT from 'transform' !!!
  scale = scene->getBound().radius() * fScaleDrag2;
  dragger2->setMatrix(osg::Matrix::scale(scale, scale, scale) *
                      osg::Matrix::translate(scene->getBound().center()));
  dragger2->addTransformUpdating(transform);
  dragger2->setHandleEvents(true);
  dragger2->setActivationKeyEvent(keyDrag2);
  // dragger2->setDraggerActive(false); // not working ?
  dragger2->setNodeMask(m_ShowRotateDragger ? 0xffffffff : 0x0);
  v_draggers2.push_back(dragger2);

  return root;
}
// osgManipulator


template< typename HalfedgeGraph, typename PointMap >
osg::ref_ptr< osg::Group >
FEVV::SimpleViewer::createMesh(
    HalfedgeGraph *_g,
    PMapsContainer *_pmaps,
    PointMap *_pm,
    std::string _mesh_file,
    osg::ref_ptr< osg::Group > _group)
{
  std::vector< osg::ref_ptr< osg::Geometry > > l_geometries, l_geometriesL,
      l_geometriesP, l_geometries_edges, l_geometries_vertices, l_geometries_normals;
  std::vector< osg::ref_ptr< osg::Vec3Array > > l_vertexArrays,
      l_vertexArrays_edges, l_vertexArrays_vertices, l_vertexArrays_normals;
  std::vector< osg::ref_ptr< osg::Vec3Array > > l_normalsArrays,
      l_normalsArraysF, l_normalsArrays_edges, l_normalsArrays_vertices, l_tangentsArrays;
  std::vector< osg::ref_ptr< osg::Vec4Array > > l_colorsArrays,
      l_colorsArrays_edges, l_colorsArrays_vertices, l_colorsArrays_normals;
  std::vector< osg::ref_ptr< osg::Vec2Array > > l_texcoordsArrays;

  osg::Geode *geode = internal_createMesh(_g,
                                          _pmaps,
                                          l_geometries,
                                          l_geometriesL,
                                          l_geometriesP,
                                          l_geometries_edges,
                                          l_geometries_vertices,
                                          l_geometries_normals,
                                          l_vertexArrays,
                                          l_vertexArrays_edges,
                                          l_vertexArrays_vertices,
                                          l_vertexArrays_normals,
                                          l_normalsArrays,
                                          l_normalsArraysF,
                                          l_normalsArrays_edges,
                                          l_normalsArrays_vertices,
                                          l_tangentsArrays,
                                          l_colorsArrays,
                                          l_colorsArrays_edges,
                                          l_colorsArrays_vertices,
                                          l_colorsArrays_normals,
                                          l_texcoordsArrays,
                                          _pm,
                                          _mesh_file);

#ifdef MANIPULATOR
  _group->addChild(addDraggersToScene(
      geode, "TabBoxDragger", 2.0, 't', "TrackballDragger", 1.0, 'r'));
#else
  _group->addChild(geode); //_group->setName("_group");
#endif

  v_geometries.push_back(l_geometries);
  v_geometriesL.push_back(l_geometriesL);
  v_geometriesP.push_back(l_geometriesP);
  v_geometries_edges.push_back(l_geometries_edges);
  v_geometries_vertices.push_back(l_geometries_vertices);
  v_geometries_normals.push_back(l_geometries_normals);
  v_vertexArrays.push_back(l_vertexArrays);
  v_vertexArrays_edges.push_back(l_vertexArrays_edges);
  v_vertexArrays_vertices.push_back(l_vertexArrays_vertices);
  v_vertexArrays_normals.push_back(l_vertexArrays_normals);
  v_normalsArrays.push_back(l_normalsArrays);
  v_normalsArraysF.push_back(l_normalsArraysF);
  v_normalsArrays_edges.push_back(l_normalsArrays_edges);
  v_normalsArrays_vertices.push_back(l_normalsArrays_vertices);
  v_tangentsArrays.push_back(l_tangentsArrays);
  v_colorsArrays.push_back(l_colorsArrays);
  v_colorsArrays_edges.push_back(l_colorsArrays_edges);
  v_colorsArrays_vertices.push_back(l_colorsArrays_vertices);
  v_colorsArrays_normals.push_back(l_colorsArrays_normals);
  v_texcoordsArrays.push_back(l_texcoordsArrays);

  v_mixed_meshes.push_back(_g);
  v_meshes_names.push_back(_mesh_file);
  v_properties_maps.push_back(_pmaps);
  v_geodes.push_back(geode);
  v_meshIsSelected.push_back(false);

  // select newly created mesh
  setNodeSelected(geode, true);

  // time
  i_time++;

  /*if (m_space_time)
  {
          if (m_time)
          {
                  current_i_time = i_time;

                  //setNodeSelected(v_geodes[current_i_time], true); // UTILE ?
          }
  }*/
  // time

  return _group;
}

template< typename HalfedgeGraph, typename PointMap >
void
FEVV::SimpleViewer::drawMesh(HalfedgeGraph *_g,
                                              PMapsContainer *_pmaps,
                                              PointMap *_pm,
                                              std::string _mesh_file)
{
  QApplication::setOverrideCursor(Qt::BusyCursor);
  SimpleWindow *sw = static_cast< SimpleWindow * >(
      getWindow()); // here static_cast instead of dynamic_cast ONLY to be
                    // consistent with 2 others for statusBar() (only for OSX
                    // and because of plugins... don't understand why...)
  sw->statusBar()->showMessage(
      QObject::tr("Draw : calculate normals & tangents...") /*, 2000*/);

  // ---
  m_redraw = false;
  m_recomputeNT_if_redraw = false;
  m_recreateOSGobj_if_redraw = true;
  // ---

  addGroup(createMesh(_g, _pmaps, _pm, _mesh_file));

  sw->statusBar()->showMessage(QObject::tr("") /*, 2000*/);
  QApplication::restoreOverrideCursor();
}

template< typename HalfedgeGraph, typename PointMap >
void
FEVV::SimpleViewer::redrawMesh(HalfedgeGraph *_g,
                                                PMapsContainer *_pmaps,
                                                PointMap *_pm,
                                                std::string _mesh_file)
{
  QApplication::setOverrideCursor(Qt::BusyCursor);
  SimpleWindow *sw = static_cast< SimpleWindow * >(
      getWindow()); // here static_cast instead of dynamic_cast only for OSX and
                    // because of plugins... don't understand why...
  if(m_recomputeNT_if_redraw)
    sw->statusBar()->showMessage(
        QObject::tr("Redraw : recalculate normals & tangents...") /*, 2000*/);
  else
    sw->statusBar()->showMessage(QObject::tr(
        "Redraw : get already calculated normals & tangents...") /*, 2000*/);

  // std::cout << "redrawMesh redrawMesh redrawMesh" << std::endl;

  // ---
  m_redraw = true;
  m_step = 0.;
  // ---

  unsigned int position;
  position = getMeshId(static_cast< void * >(_g));

  if(! Assert::check(position != -1,
        "mesh was not found. Leaving...",
        "SimpleViewer::redrawMesh"))
  {
    return;
  }

  // update mesh name during redraw if != ""
  if(_mesh_file != std::string(""))
    v_meshes_names[position] = _mesh_file;

  // MT - IMPORTANT : remove previous MATERIAL
  // v_geodes[position]->getOrCreateStateSet()->removeAttribute(osg::StateAttribute::MATERIAL);

  // MT - IMPORTANT and BETTER : reset StateSet
  v_geodes[position]->setStateSet(NULL);

  if(m_recreateOSGobj_if_redraw) // mesh
  {
    v_geometries[position].clear();
    v_geometriesL[position].clear();
    v_geometriesP[position].clear();
    v_vertexArrays[position].clear();
    v_colorsArrays[position].clear();

    v_normalsArrays[position].clear();
    v_normalsArraysF[position].clear();
    v_tangentsArrays[position].clear();

    v_texcoordsArrays[position].clear();
  }
  if(m_recreateOSGobj_if_redraw) // superimpose and 'only_pts' mode only
  {
    v_geometries_edges[position].clear();
    v_geometries_vertices[position].clear();
    v_geometries_normals[position].clear();
    v_vertexArrays_edges[position].clear();
    v_vertexArrays_vertices[position].clear();
    v_vertexArrays_normals[position].clear();
    v_colorsArrays_edges[position].clear();
    v_colorsArrays_vertices[position].clear();
    v_colorsArrays_normals[position].clear();

    v_normalsArrays_edges[position].clear();
    v_normalsArrays_vertices[position].clear();
  }

  internal_createMesh(v_geodes[position],
                      _g,
                      _pmaps,
                      v_geometries[position],
                      v_geometriesL[position],
                      v_geometriesP[position],
                      v_geometries_edges[position],
                      v_geometries_vertices[position],
                      v_geometries_normals[position],
                      v_vertexArrays[position],
                      v_vertexArrays_edges[position],
                      v_vertexArrays_vertices[position],
                      v_vertexArrays_normals[position],
                      v_normalsArrays[position],
                      v_normalsArraysF[position],
                      v_normalsArrays_edges[position],
                      v_normalsArrays_vertices[position],
                      v_tangentsArrays[position],
                      v_colorsArrays[position],
                      v_colorsArrays_edges[position],
                      v_colorsArrays_vertices[position],
                      v_colorsArrays_normals[position],
                      v_texcoordsArrays[position],
                      _pm,
                      v_meshes_names[position]);

  sw->statusBar()->showMessage(QObject::tr("") /*, 2000*/);
  QApplication::restoreOverrideCursor();
}

template< typename HalfedgeGraph >
void
FEVV::SimpleViewer::centerMesh(HalfedgeGraph *_g)
{
  unsigned int position;
  position = getMeshId(static_cast< void * >(_g));

  if(! Assert::check(position != -1,
        "mesh was not found. Leaving...",
        "SimpleViewer::centerMesh"))
  {
    return;
  }

  // ---

  // std::cout << "getNumViews() : " << getNumViews() << std::endl;
  // osgViewer::View* _osgView = getView(0); // for osgViewer::CompositeViewer
  osgViewer::View *_osgView =
      dynamic_cast< osgViewer::View * >(this); // for osgViewer::Viewer
  // osgViewer::View* _osgView = getViewWithFocus();

  if(_osgView)
  {
#ifdef MANIPULATOR
    _osgView->getCameraManipulator()->setNode(v_geodes[position]->getParent(0));


    // ELO-note: v_geodes[position]->getParent(0) is an osg::MatrixTransform
    //           which is a group with an osg::Matrix
    //           see
    //           http://camlosg.sourceforge.net/osg/classosg_1_1MatrixTransform.html


    // std::cout << "centerMesh (MANIPULATOR) \"" <<
    // v_geodes[position]/*->getParent(0)*/->getName() << "\"" << std::endl;
#else
    _osgView->getCameraManipulator()->setNode(v_geodes[position]);
    // std::cout << "centerMesh \"" << v_geodes[position]->getName() << "\"" <<
    // std::endl;
#endif
    _osgView->getCameraManipulator()->computeHomePosition();
    _osgView->home();


    // FEVV::Debug::print_osg_tree_from_node(v_geodes[position]->getParent(0)->getParent(0));
  }
}


inline
osg::Matrix
FEVV::SimpleViewer::getTransformMatrixOsg(
    unsigned int position)
{
  assert(position < v_geodes.size());
  osg::MatrixTransform *grp_MatrixTransform =
      dynamic_cast< osg::MatrixTransform * >(v_geodes[position]->getParent(0));
  assert(grp_MatrixTransform != nullptr);
  osg::Matrix matrix = grp_MatrixTransform->getMatrix();

  return matrix; // 4x4 homogeneous matrix
}


inline
Eigen::Matrix4d
FEVV::SimpleViewer::getTransformMatrixEigen(
    unsigned int position)
{
  osg::Matrix osg_mat = getTransformMatrixOsg(position);

  // convert OSG transform matrix to Eigen matrix
  // transposition needed!
  Eigen::Matrix4d eigen_mat;
  eigen_mat << osg_mat(0, 0), osg_mat(1, 0), osg_mat(2, 0), osg_mat(3, 0),
      osg_mat(0, 1), osg_mat(1, 1), osg_mat(2, 1), osg_mat(3, 1), osg_mat(0, 2),
      osg_mat(1, 2), osg_mat(2, 2), osg_mat(3, 2), osg_mat(0, 3), osg_mat(1, 3),
      osg_mat(2, 3), osg_mat(3, 3);

  // std::cout << "eigen_mat = \n" << eigen_mat << std::endl;

  return eigen_mat; // 4x4 homogeneous matrix
}


inline
void
FEVV::SimpleViewer::resetTransformMatrix(unsigned int position)
{
  assert(position < v_geodes.size());
  osg::MatrixTransform *grp_MatrixTransform =
      dynamic_cast< osg::MatrixTransform * >(v_geodes[position]->getParent(0));
  assert(grp_MatrixTransform != nullptr);

  osg::Matrix identity;
  identity.makeIdentity();
  grp_MatrixTransform->setMatrix(identity);
}


template< typename HalfedgeGraph >
void
FEVV::SimpleViewer::draw_or_redraw_mesh(
    /*const */ HalfedgeGraph *_g,
    /*const */ PMapsContainer *_pmaps,
    bool _redraw,
    bool _recomputeNT_if_redraw,
    std::string _mesh_filename,
    bool _recreateOSGobj_if_redraw,
    float _step)
{
  auto pm = get(boost::vertex_point, *_g);

  if(!_redraw)
  {
    m_step = _step;

    // draw mesh
    drawMesh(_g,
             _pmaps,
             &pm,           /*point map*/
             _mesh_filename /*mesh filename*/
    );

    centerMesh(_g);
  }
  else
  {
    m_recomputeNT_if_redraw = _recomputeNT_if_redraw;
    m_recreateOSGobj_if_redraw = _recreateOSGobj_if_redraw;

    // redraw mesh
    redrawMesh(_g,
               _pmaps,
               &pm,           /*point map*/
               _mesh_filename /*mesh filename*/
    );
  }
}


inline
void
FEVV::SimpleViewer::activate_time_mode()
{
  SimpleWindow *sw = static_cast< SimpleWindow * >(
      getWindow()); // here static_cast instead of dynamic_cast only for OSX and
                    // because of plugins... don't understand why...

  sw->activate_time_mode();
}


inline
void
FEVV::SimpleViewer::activate_space_mode()
{
  SimpleWindow *sw = static_cast< SimpleWindow * >(
      getWindow()); // here static_cast instead of dynamic_cast only for OSX and
                    // because of plugins... don't understand why...

  sw->activate_space_mode();
}


inline
void
FEVV::SimpleViewer::updateSWModelList()
{
  SimpleWindow *sw = static_cast< SimpleWindow * >(
      getWindow()); // here static_cast instead of dynamic_cast only for OSX and
                    // because of plugins... don't understand why...

  sw->update();
}


inline
void
FEVV::SimpleViewer::setNodeSelected(osg::Node *_geode,
                                                     bool _isSelected)
{
  osg::Geode *geode = dynamic_cast< osg::Geode * >(_geode);
  if(geode != nullptr)
  {
    unsigned int position = 0;
    for(osg::Geode *g : v_geodes)
    {
      if(g == geode)
      {
        v_meshIsSelected[position] = _isSelected;

        if(m_space_time)
        {
          if(_isSelected)
          {
            current_i_time = position; // IMPORTANT !!!

            if(m_time)
            {
              for(unsigned i = 0; i < v_geodes.size(); i++)
                v_geodes[i]->setNodeMask(0x0);
              for(unsigned i = 0; i < v_draggers1.size(); i++)
                v_draggers1[i]->setNodeMask(0x0);
              for(unsigned i = 0; i < v_draggers2.size(); i++)
                v_draggers2[i]->setNodeMask(0x0);

              v_geodes[position]->setNodeMask(0xffffffff);
              if(v_draggers1.size())
                v_draggers1[position]->setNodeMask(
                    m_ShowTranslateDragger ? 0xffffffff : 0x0);
              if(v_draggers2.size())
                v_draggers2[position]->setNodeMask(
                    m_ShowRotateDragger ? 0xffffffff : 0x0);
            }
          }
        }

        break;
      }
      ++position;
    }
  }

#if 0
    osg::ref_ptr< osgFX::Outline > outline = new osgFX::Outline;
    if ( _isSelected )
    {
        outline->setWidth(3);
        outline->setColor(Helpers::ColorConverter(Color::Yellow()));
        outline->addChild(_geode);
        for( auto parent : _geode->getParents() )
        {
            if( parent == outline )
            {
                continue;
            }
            parent->replaceChild( _geode, outline );
        }

        return;
    }
    else
    {
        for( auto parent : _geode->getParents() )
        {
            if( parent->isSameKindAs( outline ))
            {
                for( auto grandpa : parent->getParents() )
                {
                    grandpa->replaceChild( parent, _geode );
                }
                break;
            }
        }
    }
    return;
#endif

#if 0
  //TODO-elo-DEBUG
  std::cout << "viewer " << this << "  v_meshIsSelected =";
  for(bool b: v_meshIsSelected)
    std::cout << " " << b;
  std::cout << std::endl;
#endif
}


inline
bool
FEVV::SimpleViewer::isNodeSelected(osg::Node *_geode)
{
  osg::Geode *geode = dynamic_cast< osg::Geode * >(_geode);
  if(geode != nullptr)
  {
    unsigned int position = 0;
    for(osg::Geode *g : v_geodes)
    {
      if(g == geode)
      {
        return v_meshIsSelected[position];
      }
      ++position;
    }
  }

  return false;
}


inline
size_t
FEVV::SimpleViewer::getMeshId(const void *mesh_ptr)
{
  {
    for(size_t i_pos = 0; i_pos < v_mixed_meshes.size(); i_pos++)
    {
      if(v_mixed_meshes[i_pos].first == mesh_ptr)
        return i_pos;
    }

    return -1; // not found
  }
}
