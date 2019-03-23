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

template< typename HalfedgeGraph >
FEVV::SimpleViewer< HalfedgeGraph >::SimpleViewer() : BaseViewerOSG()
{
}

template< typename HalfedgeGraph >
FEVV::SimpleViewer< HalfedgeGraph >::~SimpleViewer()
{
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

      SimpleWindow *sw = dynamic_cast< SimpleWindow * >(getWindow());
      sw->update();
    }
  }

  /*for (unsigned int ii = 0; ii < v_geodes.size(); ++ii)
          delete v_geodes[ii]; // -> impossible, car "smart pointer"*/

  /*for (unsigned int ii = 0; ii < v_draggers1.size(); ++ii)
          delete v_draggers1[ii];	// -> impossible, car "smart pointer"
  for (unsigned int ii = 0; ii < v_draggers2.size(); ++ii)
          delete v_draggers2[ii];	// -> impossible, car "smart pointer"*/

  for(unsigned int ii = 0; ii < v_meshes.size(); ++ii)
  {
    {
      delete v_meshes[ii];
      // std::cout << "--> delete v_meshes[ii] (" <<
      // bavQt->windowTitle().toStdString() << ") in ~SimpleViewer" << std::endl;
    }
  }

  for(unsigned int ii = 0; ii < v_properties_maps.size(); ++ii)
  {
    delete v_properties_maps[ii];
    // std::cout << "--> delete v_properties_maps[ii] (" <<
    // bavQt->windowTitle().toStdString() << ") in ~SimpleViewer" << std::endl;
  }
}

template< typename HalfedgeGraph >
void
FEVV::SimpleViewer< HalfedgeGraph >::init()
{
  if(!Assert::check(
         !bIsInit, "is already init. Leaving...", "SimpleViewer::init"))
  {
    return;
  }

  // disable the default setting of viewer.done() by pressing Escape.
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
  addGroup(gizmo);
  grid = Debug::createUnitGrid();
  addGroup(grid);

  // hud = Debug::createHud(hudText);
  // addGroup(hud);

  // time
  i_time = current_i_time = -1;
  // time

  std::cout << "SimpleViewer init" << std::endl;

  bIsInit = true;
}

template< typename HalfedgeGraph >
bool
FEVV::SimpleViewer< HalfedgeGraph >::isInit() const
{
  return bIsInit;
}

template< typename HalfedgeGraph >
bool
FEVV::SimpleViewer< HalfedgeGraph >::isValid() const
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

template< typename HalfedgeGraph >
bool
FEVV::SimpleViewer< HalfedgeGraph >::changeBackgroundColor(
    const FEVV::Color &_color)
{
  // osgViewer::View* _osgView = getView(0); // for osgViewer::CompositeViewer
  osgViewer::View *_osgView =
      dynamic_cast< osgViewer::View * >(this); // for osgViewer::Viewer
  _osgView->getCamera()->setClearColor(Helpers::ColorConverter(_color));

  return true;
}

template< typename HalfedgeGraph >
bool
FEVV::SimpleViewer< HalfedgeGraph >::saveScreenshot(const std::string &_name)
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

template< typename HalfedgeGraph >
void
FEVV::SimpleViewer< HalfedgeGraph >::addModel(Model *_geode)
{
  root_node->addChild(_geode);
  if(myWindow != nullptr)
  {
    myWindow->notify();
  }

  // osgUtil::Optimizer optimizer;
  // optimizer.optimize( root_node );
}

template< typename HalfedgeGraph >
void
FEVV::SimpleViewer< HalfedgeGraph >::addGroup(Group *_group)
{
  root_node->addChild(_group);
  if(myWindow != nullptr)
  {
    myWindow->notify();
  }

  // osgUtil::Optimizer optimizer;
  // optimizer.optimize( root_node );
}

template< typename HalfedgeGraph >
typename FEVV::SimpleViewer< HalfedgeGraph >::DataModelVector *
FEVV::SimpleViewer< HalfedgeGraph >::getDataModel()
{
  visitor->reset();

  root_node->traverse(*visitor);

  return visitor->exportResults();
}

template< typename HalfedgeGraph >
std::vector< osg::Geode * >
FEVV::SimpleViewer< HalfedgeGraph >::getSelectedGeodes()
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

template< typename HalfedgeGraph >
std::vector< osg::Geode * >
FEVV::SimpleViewer< HalfedgeGraph >::getGeodes()
{
  return v_geodes;
}

template< typename HalfedgeGraph >
std::vector< HalfedgeGraph * >
FEVV::SimpleViewer< HalfedgeGraph >::getSelectedMeshes()
{
  std::vector< HalfedgeGraph * > result;
  unsigned int i_pos = 0;
  for(bool b : v_meshIsSelected)
  {
    if(b)
    {
      result.push_back(v_meshes[i_pos]);
    }
    ++i_pos;
  }

  return result;
}

template< typename HalfedgeGraph >
std::vector< HalfedgeGraph * >
FEVV::SimpleViewer< HalfedgeGraph >::getMeshes()
{
  return v_meshes;
}

template< typename HalfedgeGraph >
std::vector< std::string >
FEVV::SimpleViewer< HalfedgeGraph >::getSelectedMeshesNames()
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

template< typename HalfedgeGraph >
std::vector< std::string >
FEVV::SimpleViewer< HalfedgeGraph >::getMeshesNames()
{
  return v_meshes_names;
}

template< typename HalfedgeGraph >
std::vector< FEVV::PMapsContainer * >
FEVV::SimpleViewer< HalfedgeGraph >::getSelected_properties_maps()
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

template< typename HalfedgeGraph >
std::vector< FEVV::PMapsContainer * >
FEVV::SimpleViewer< HalfedgeGraph >::get_properties_maps()
{
  return v_properties_maps;
}

template< typename HalfedgeGraph >
std::vector< osg::Group * >
FEVV::SimpleViewer< HalfedgeGraph >::getSelectedDraggers1()
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

template< typename HalfedgeGraph >
std::vector< osg::Group * >
FEVV::SimpleViewer< HalfedgeGraph >::getDraggers1()
{
  return v_draggers1;
}

template< typename HalfedgeGraph >
std::vector< osg::Group * >
FEVV::SimpleViewer< HalfedgeGraph >::getSelectedDraggers2()
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

template< typename HalfedgeGraph >
std::vector< osg::Group * >
FEVV::SimpleViewer< HalfedgeGraph >::getDraggers2()
{
  return v_draggers2;
}

template< typename HalfedgeGraph >
HalfedgeGraph *
FEVV::SimpleViewer< HalfedgeGraph >::getMesh(unsigned int _position)
{
  Assert::check(_position < v_meshes.size(),
                "The given position must be lower than v_meshes.size()",
                "SimpleViewer::getMesh(int)");

  return v_meshes[_position];
}

template< typename HalfedgeGraph >
template< typename PointMap,
          typename VertexNormalMap,
          typename FaceNormalMap,
          typename VertexColorMap,
          typename FaceColorMap,
          typename VertexUVMap,
          typename HalfedgeUVMap >
inline osg::Geode *
FEVV::SimpleViewer< HalfedgeGraph >::internal_createMesh(
    HalfedgeGraph *_g,
    PMapsContainer *_pmaps,
    std::map< vertex_descriptor, unsigned int > &_mapVertex,
    std::map< face_descriptor, unsigned int > &_mapFace,
    PointMap *_pm,
    VertexNormalMap *_vt_nm,
    FaceNormalMap *_f_nm,
    VertexColorMap *_vt_cm,
    FaceColorMap *_f_cm,
    VertexUVMap *_vt_uv_m,
    HalfedgeUVMap *_het_uv_m,
    int _texture_type,
    std::string _texture_file,
    std::string _mesh_file)
{
  osg::Geode *geode =
      new osg::Geode; // MT OSG_GEODE : bug from JL ? Ici devrait être
                      // osg::ref_ptr<osg::Geode> geode = new osg::Geode(); ???
  internal_createMesh(geode,
                      _g,
                      _pmaps,
                      _mapVertex,
                      _mapFace,
                      _pm,
                      _vt_nm,
                      _f_nm,
                      _vt_cm,
                      _f_cm,
                      _vt_uv_m,
                      _het_uv_m,
                      _texture_type,
                      _texture_file,
                      _mesh_file);
  return geode;
}

template< typename HalfedgeGraph >
template< typename PointMap,
          typename VertexNormalMap,
          typename FaceNormalMap,
          typename VertexColorMap,
          typename FaceColorMap,
          typename VertexUVMap,
          typename HalfedgeUVMap >
void
FEVV::SimpleViewer< HalfedgeGraph >::internal_createMesh(
    osg::Geode *&geode,
    HalfedgeGraph *_g,
    PMapsContainer *_pmaps,
    std::map< vertex_descriptor, unsigned int > &_mapVertex,
    std::map< face_descriptor, unsigned int > &_mapFace,
    PointMap *_pm,
    VertexNormalMap *_vt_nm,
    FaceNormalMap *_f_nm,
    VertexColorMap *_vt_cm,
    FaceColorMap *_f_cm,
    VertexUVMap *_vt_uv_m,
    HalfedgeUVMap *_het_uv_m,
    int _texture_type,
    std::string _texture_file,
    std::string _mesh_file)
{
  // property maps stuff
  using FaceMaterialMap =
      typename FEVV::PMap_traits< FEVV::face_material_t,
                                  HalfedgeGraph >::pmap_type;
  using MeshMaterialsMap =
      typename FEVV::PMap_traits< FEVV::mesh_materials_t,
                                  HalfedgeGraph >::pmap_type;
  using MeshGuipropertiesMap =
      typename FEVV::PMap_traits< FEVV::mesh_guiproperties_t,
                                  HalfedgeGraph >::pmap_type;
  // RM: tangents are used primarily for normal mapping
  using VertexTangentMap =
      typename FEVV::PMap_traits< FEVV::vertex_tangent_t,
                                  HalfedgeGraph >::pmap_type;

  VertexNormalMap v_nm;
  FaceNormalMap f_nm;
  VertexColorMap v_cm;
  FaceColorMap f_cm;
  VertexUVMap v_uvm;
  HalfedgeUVMap h_uvm;

  VertexTangentMap vt_m;
  FaceMaterialMap f_mm;
  MeshMaterialsMap m_mm;
  MeshGuipropertiesMap m_gpm;

  _vt_nm = nullptr;
  _f_nm = nullptr;
  _vt_cm = nullptr;
  _f_cm = nullptr;
  _vt_uv_m = nullptr;
  _het_uv_m = nullptr;

  VertexTangentMap *v_tan_m = nullptr;
  FaceMaterialMap *_f_mm = nullptr;
  MeshMaterialsMap *_m_mm = nullptr;
  MeshGuipropertiesMap *_m_gpm = nullptr;
  size_t _m_mm_size = 0;

  // textures stuff
  _texture_type = NO_TEXCOORDS;
  _texture_file = "";
  // TODO-elo  remove the '_texture_type' variable ;
  //          the texture type can be deduced of the property maps
  //          that are read from the file ;

  FEVV::Filters::translate(
      *_g, *_pm, m_step, 0., 0.); // TEMP for test

  // retrieve or create mesh gui-properties property map
  if(has_map(*_pmaps, FEVV::mesh_guiproperties))
  {
    m_gpm = get_property_map(FEVV::mesh_guiproperties, *_g, *_pmaps);
    _m_gpm = &m_gpm;
    //std::cout << "[SimpleViewer] **************get mesh_guiproperties property_map" << std::endl;
  }
  else
  {
    m_gpm = make_property_map(FEVV::mesh_guiproperties, *_g);
    FEVV::Types::GuiProperties gui_props;
    put(m_gpm, 0, gui_props);
    put_property_map(FEVV::mesh_guiproperties, *_g, *_pmaps, m_gpm);
    _m_gpm = &m_gpm;
    //std::cout << "[SimpleViewer] **************make mesh_guiproperties property_map" << std::endl;
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

  if(m_redraw) // ONLY IF REDRAW via GUI or CODE
  {
   if( (m_UseVertexColor) || (m_UseFaceColor) || (m_UseTexture) ) // else automatic detection, as for a first DRAW
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
      //_texture_file = "";

      _f_mm = nullptr;
      _m_mm = nullptr;
      _m_mm_size = 0;
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
                // -- TODO MT - TEMP - HERE WE FORCE a full RED vertex-color map
        for the REDRAW using Vector = typename GeometryTraits::Vector;

                auto iterator_pair = vertices(*_g); // vertices() returns a
        vertex_iterator pair vertex_iterator vi = iterator_pair.first;
                vertex_iterator vi_end = iterator_pair.second;
                for (; vi != vi_end; ++vi)
                {
                        v_cm[*vi] = Vector(1., 0., 0.); // RGB
                }

                _vt_cm = &v_cm;
                // -- TODO MT - TEMP - HERE WE FORCE a full RED vertex-color map
        for the REDRAW
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
            _texture_type = NO_TEXCOORDS; // if face material map is missing, we
                                          // can not display any texture in
                                          // HALFEDGE_TEXCOORDS2D mode
        }

        if(has_map(*_pmaps, FEVV::mesh_materials))
        {
          m_mm = get_property_map(FEVV::mesh_materials, *_g, *_pmaps);
          _m_mm = &m_mm;

          _m_mm_size = std::distance(m_mm.storage_begin(), m_mm.storage_end());
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

  _mapVertex.clear();
  _mapFace.clear();

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
  std::vector< osg::ref_ptr< osg::Geometry > > geometries, geometries_edges,
      geometries_vertices;
  std::vector< osg::ref_ptr< osg::Vec3Array > > vertexArrays,
      vertexArrays_edges, vertexArrays_vertices;
  std::vector< osg::ref_ptr< osg::Vec3Array > > normalsArrays;
  std::vector< osg::ref_ptr< osg::Vec3Array > > tangentsArrays;
  std::vector< osg::ref_ptr< osg::Vec4Array > > colorsArrays,
      colorsArrays_edges, colorsArrays_vertices;
  std::vector< osg::ref_ptr< osg::Vec2Array > > texcoordsArrays;

  geometries.push_back(new osg::Geometry()); geometries[0]->setSupportsDisplayList(true); geometries[0]->setUseDisplayList(true); geometries[0]->setUseVertexBufferObjects(false); /*geometries[0]->setUseVertexArrayObject(false);*/
  geometries_edges.push_back(new osg::Geometry()); geometries_edges[0]->setSupportsDisplayList(true); geometries_edges[0]->setUseDisplayList(true); geometries_edges[0]->setUseVertexBufferObjects(false); /*geometries_edges[0]->setUseVertexArrayObject(false);*/
  geometries_vertices.push_back(new osg::Geometry()); geometries_vertices[0]->setSupportsDisplayList(true); geometries_vertices[0]->setUseDisplayList(true); geometries_vertices[0]->setUseVertexBufferObjects(false); /*geometries_vertices[0]->setUseVertexArrayObject(false);*/
  vertexArrays.push_back(new osg::Vec3Array);
  vertexArrays_edges.push_back(new osg::Vec3Array);
  vertexArrays_vertices.push_back(new osg::Vec3Array);
  normalsArrays.push_back(new osg::Vec3Array);
  tangentsArrays.push_back(new osg::Vec3Array);
  colorsArrays.push_back(new osg::Vec4Array);
  colorsArrays_edges.push_back(new osg::Vec4Array);
  colorsArrays_vertices.push_back(new osg::Vec4Array);
  texcoordsArrays.push_back(new osg::Vec2Array);
  // std::cout << "---------> geometries.push_back" << std::endl;
  size_t mtl_id = 0;

  // add ONLY if multi-textures
  for(size_t mi = 1; mi < _m_mm_size; mi++)
  {
    geometries.push_back(new osg::Geometry()); geometries[mi]->setSupportsDisplayList(true); geometries[mi]->setUseDisplayList(true); geometries[mi]->setUseVertexBufferObjects(false); /*geometries[mi]->setUseVertexArrayObject(false);*/
    geometries_edges.push_back(new osg::Geometry()); geometries_edges[mi]->setSupportsDisplayList(true); geometries_edges[mi]->setUseDisplayList(true); geometries_edges[mi]->setUseVertexBufferObjects(false); /*geometries_edges[mi]->setUseVertexArrayObject(false);*/
    geometries_vertices.push_back(new osg::Geometry()); geometries_vertices[mi]->setSupportsDisplayList(true); geometries_vertices[mi]->setUseDisplayList(true); geometries_vertices[mi]->setUseVertexBufferObjects(false); /*geometries_vertices[mi]->setUseVertexArrayObject(false);*/
    vertexArrays.push_back(new osg::Vec3Array);
    vertexArrays_edges.push_back(new osg::Vec3Array);
    vertexArrays_vertices.push_back(new osg::Vec3Array);
    normalsArrays.push_back(new osg::Vec3Array);
    tangentsArrays.push_back(new osg::Vec3Array);
    colorsArrays.push_back(new osg::Vec4Array);
    colorsArrays_edges.push_back(new osg::Vec4Array);
    colorsArrays_vertices.push_back(new osg::Vec4Array);
    texcoordsArrays.push_back(new osg::Vec2Array);
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

  //std::cout << "---------> texture_corner_mode_on: " << texture_corner_mode_on << std::endl;
  //std::cout << "---------> texture_vertex_mode_on: " << texture_vertex_mode_on << std::endl;

  /// Adding edges - superimpose only
  if(m_RenderSuperimposedEdges)
  {
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
      vertexArrays_edges[mtl_id]->push_back(
          Helpers::VectorConverter< HalfedgeGraph >(ps));

      // retrieve vertex (target)
      vt = target(*eb, *_g);

      // retrieve vertex-point (target)
      pt = (*_pm)[vt];
      vertexArrays_edges[mtl_id]->push_back(
          Helpers::VectorConverter< HalfedgeGraph >(pt));

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
    }

    geometries_edges[mtl_id]->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINES, 0, vertexArrays_edges[mtl_id]->size()));

	// set line width
	osg::ref_ptr< osg::LineWidth > linewidth = new osg::LineWidth();
	linewidth->setWidth(3.0f);
	geometries_edges[mtl_id]
	  ->getOrCreateStateSet()
	  ->setAttribute /*setAttributeAndModes*/ (linewidth,
	                                           osg::StateAttribute::ON);

	// light
	geometries_edges[mtl_id]->getOrCreateStateSet()->setMode(
	  GL_LIGHTING,
	  osg::StateAttribute::OFF); // light always OFF for superimpose edges
  }

  /// Adding vertices - superimpose and 'only_pts' mode only
  size_t nb_faces = size_of_faces(*_g);
  if(m_RenderSuperimposedVertices || m_RenderSuperimposedVertices_Big || (nb_faces==0)) // last test for 'only_pts' mode
  {
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

      sizeSPoints++;

      // color
      if(_vt_cm)
        colorsArrays_vertices[mtl_id]->push_back(
            Helpers::VectorColorConverter< HalfedgeGraph >(
                get(v_cm, *v_it))); // user/filter/plugin colors
      else
        colorsArrays_vertices[mtl_id]->push_back(
            Helpers::ColorConverter(Color::Green())); // default color
    }

    geometries_vertices[mtl_id]->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POINTS, 0, vertexArrays_vertices[mtl_id]->size()));

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

    // ---

    _vt_cm = SAVE_vt_cm;
  }

  /// Adding faces
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
        if(m_SmoothFlat_Shading)
        {
          normal = (*_vt_nm)[vd.back()];
          normalsArrays[mtl_id]->push_back(
              Helpers::VectorConverter< HalfedgeGraph >(
                  normal)); // smooth mode
        }
        else
        {
          normal = (*_f_nm)[*fb];
          normalsArrays[mtl_id]->push_back(
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
        //std::cout << "---> QUAD-QUAD-QUAD-QUAD-QUAD" << std::endl;
      }
      else
      {
        //std::cout << "---> POLY-POLY-POLY-POLY-POLY" << std::endl;
      }
    }

    geometries[mtl_id]->addPrimitiveSet(new osg::DrawArrays(drawing_method, vertexArrays[mtl_id]->size() - num_vertices_in_face, num_vertices_in_face));

    // populate face normal array
    if(_vt_nm != nullptr) // normal per vertex (see above) -> already calculated
    {
      // do nothing
    }
    else if(_f_nm != nullptr) // normal per face -> already calculated
    {
      normal = (*_f_nm)[*fb];
      normalsArrays[mtl_id]->push_back(
          Helpers::VectorConverter< HalfedgeGraph >(normal));
    }
    else // if none - BUT NEVER HAPPENS - re-compute normal for current face
    {
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
    }

    // populate face color array
    if(_f_cm != nullptr)
    {
      colorsArrays[mtl_id]->push_back(
          Helpers::VectorColorConverter< HalfedgeGraph >((*_f_cm)[*fb]));
    }

    ++sizeFace;
  }

  sw->statusBar()->showMessage(QObject::tr("") /*, 2000*/);
  QApplication::restoreOverrideCursor();

  std::cout << "[SimpleViewer] I have drawn " << sizeVertex << " vertices and "
            << sizeFace << " faces." << std::endl;
  std::cout << "[SimpleViewer] I have also drawn " << sizeSPoints << " (superimpose) points and "
            << sizeSLines << " superimpose lines." << std::endl;

  //auto gui_props = get((*_m_gpm), 0);
  //geode->setNodeMask(gui_props.is_visible ? 0xffffffff : 0x0); // 19/03/19

  const auto loadingStartTime = std::chrono::system_clock::now();

  QApplication::setOverrideCursor(Qt::BusyCursor);
  sw->statusBar()->showMessage(QObject::tr("Render mesh...") /*, 2000*/);

  if(m_RenderMode == RenderMode::RENDER_SHADERS_DIRECT_LIGHTING ||
     m_RenderMode == RenderMode::RENDER_SHADERS_INDIRECT_LIGHTING)
  {
    internal_loadShadedMesh(geode,
                            _g,
                            geometries,
                            geometries_edges,
                            geometries_vertices,
                            vertexArrays,
                            vertexArrays_edges,
                            vertexArrays_vertices,
                            normalsArrays,
                            tangentsArrays,
                            texcoordsArrays,
                            colorsArrays,
                            colorsArrays_edges,
                            colorsArrays_vertices,
                            _m_mm_size,
                            _vt_nm,
                            v_tan_m,
                            _vt_cm,
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
                            geometries,
                            geometries_edges,
                            geometries_vertices,
                            vertexArrays,
                            vertexArrays_edges,
                            vertexArrays_vertices,
                            normalsArrays,
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

  std::cout << "[SimpleViewer] Done loading mesh after "
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
  geode->setName(std::string("Mesh '") + _mesh_file + std::string("' [") +
                 bavQt->windowTitle().toStdString() + std::string("]"));
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

inline osgManipulator::Dragger *
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

template< typename HalfedgeGraph >
osg::Node *
FEVV::SimpleViewer< HalfedgeGraph >::addDraggersToScene(
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

template< typename HalfedgeGraph >
template< typename PointMap,
          typename VertexNormalMap,
          typename FaceNormalMap,
          typename VertexColorMap,
          typename FaceColorMap,
          typename VertexUVMap,
          typename HalfedgeUVMap >
osg::ref_ptr< osg::Group >
FEVV::SimpleViewer< HalfedgeGraph >::createMesh(
    HalfedgeGraph *_g,
    PMapsContainer *_pmaps,
    PointMap *_pm,
    VertexNormalMap *_vt_nm,
    FaceNormalMap *_f_nm,
    VertexColorMap *_vt_cm,
    FaceColorMap *_f_cm,
    VertexUVMap *_vt_uv_m,
    HalfedgeUVMap *_het_uv_m,
    int _texture_type,
    std::string _texture_file,
    std::string _mesh_file,
    osg::ref_ptr< osg::Group > _group)
{
  std::map< vertex_descriptor, unsigned int > mapVertex;
  std::map< face_descriptor, unsigned int > mapFace;

  osg::Geode *geode = internal_createMesh(_g,
                                          _pmaps,
                                          mapVertex,
                                          mapFace,
                                          _pm,
                                          _vt_nm,
                                          _f_nm,
                                          _vt_cm,
                                          _f_cm,
                                          _vt_uv_m,
                                          _het_uv_m,
                                          _texture_type,
                                          _texture_file,
                                          _mesh_file);

#ifdef MANIPULATOR
  _group->addChild(addDraggersToScene(
      geode, "TabBoxDragger", 2.0, 't', "TrackballDragger", 1.0, 'r'));
#else
  _group->addChild(geode); //_group->setName("_group");
#endif

  v_mapVertex.push_back(mapVertex);
  v_mapFace.push_back(mapFace);
  v_meshes.push_back(_g);
  v_meshes_names.push_back(_mesh_file);
  v_properties_maps.push_back(_pmaps);
  v_geodes.push_back(geode);
  v_meshIsSelected.push_back(false);

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

template< typename HalfedgeGraph >
template< typename PointMap,
          typename VertexNormalMap,
          typename FaceNormalMap,
          typename VertexColorMap,
          typename FaceColorMap,
          typename VertexUVMap,
          typename HalfedgeUVMap >
void
FEVV::SimpleViewer< HalfedgeGraph >::drawMesh(HalfedgeGraph *_g,
                                              PMapsContainer *_pmaps,
                                              PointMap *_pm,
                                              VertexNormalMap *_vt_nm,
                                              FaceNormalMap *_f_nm,
                                              VertexColorMap *_vt_cm,
                                              FaceColorMap *_f_cm,
                                              VertexUVMap *_vt_uv_m,
                                              HalfedgeUVMap *_het_uv_m,
                                              int _texture_type,
                                              std::string _texture_file,
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
  // ---

  addGroup(createMesh(_g,
                      _pmaps,
                      _pm,
                      _vt_nm,
                      _f_nm,
                      _vt_cm,
                      _f_cm,
                      _vt_uv_m,
                      _het_uv_m,
                      _texture_type,
                      _texture_file,
                      _mesh_file));

  sw->statusBar()->showMessage(QObject::tr("") /*, 2000*/);
  QApplication::restoreOverrideCursor();
}

template< typename HalfedgeGraph >
template< typename PointMap,
          typename VertexNormalMap,
          typename FaceNormalMap,
          typename VertexColorMap,
          typename FaceColorMap,
          typename VertexUVMap,
          typename HalfedgeUVMap >
void
FEVV::SimpleViewer< HalfedgeGraph >::redrawMesh(bool _flushMesh,
                                                HalfedgeGraph *_g,
                                                PMapsContainer *_pmaps,
                                                PointMap *_pm,
                                                VertexNormalMap *_vt_nm,
                                                FaceNormalMap *_f_nm,
                                                VertexColorMap *_vt_cm,
                                                FaceColorMap *_f_cm,
                                                VertexUVMap *_vt_uv_m,
                                                HalfedgeUVMap *_het_uv_m,
                                                int _texture_type,
                                                std::string _texture_file,
                                                std::string _mesh_file)
{
  QApplication::setOverrideCursor(Qt::BusyCursor);
  SimpleWindow *sw = static_cast< SimpleWindow * >(
      getWindow()); // here static_cast instead of dynamic_cast only for OSX and
                    // because of plugins... don't understand why...
  if (m_recomputeNT_if_redraw)
    sw->statusBar()->showMessage(
        QObject::tr("Redraw : recalculate normals & tangents...") /*, 2000*/);
  else
    sw->statusBar()->showMessage(
        QObject::tr("Redraw : get already calculated normals & tangents...") /*, 2000*/);

  // std::cout << "redrawMesh redrawMesh redrawMesh" << std::endl;

  // ---
  m_redraw = true;
  m_step = 0.;
  // ---

  unsigned int position;
  {
    unsigned int i_pos = 0;
    for(HalfedgeGraph *m : v_meshes)
    {
      if(m == _g)
      {
        position = i_pos;
        break;
      }
      ++i_pos;
    }
    if(!Assert::check(i_pos != v_meshes.size(),
                      "mesh was not found. Leaving...",
                      "SimpleViewer::redrawMesh"))
    {
      return;
    }
  }

  // update mesh name during redraw if != ""
  if(_mesh_file != std::string(""))
    v_meshes_names[position] = _mesh_file;

  if(_flushMesh) // NOW always ON so delete this param in the future
  {
    // MT - IMPORTANT : remove previous MATERIAL
    // v_geodes[position]->getOrCreateStateSet()->removeAttribute(osg::StateAttribute::MATERIAL);

    // MT - IMPORTANT and BETTER : reset StateSet
    v_geodes[position]->setStateSet(NULL);

    v_mapFace[position].clear();
    v_mapVertex[position].clear();

    internal_createMesh(v_geodes[position],
                        _g,
                        _pmaps,
                        v_mapVertex[position],
                        v_mapFace[position],
                        _pm,
                        _vt_nm,
                        _f_nm,
                        _vt_cm,
                        _f_cm,
                        _vt_uv_m,
                        _het_uv_m,
                        _texture_type,
                        _texture_file,
                        v_meshes_names[position]);
  }

  sw->statusBar()->showMessage(QObject::tr("") /*, 2000*/);
  QApplication::restoreOverrideCursor();
}

template< typename HalfedgeGraph >
void
FEVV::SimpleViewer< HalfedgeGraph >::centerMesh(HalfedgeGraph *_g)
{
  unsigned int position;
  {
    unsigned int i_pos = 0;
    for(HalfedgeGraph *m : v_meshes)
    {
      if(m == _g)
      {
        position = i_pos;
        break;
      }
      ++i_pos;
    }
    if(!Assert::check(i_pos != v_meshes.size(),
                      "mesh was not found. Leaving...",
                      "SimpleViewer::centerMesh"))
    {
      return;
    }
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
    //           see http://camlosg.sourceforge.net/osg/classosg_1_1MatrixTransform.html


    // std::cout << "centerMesh (MANIPULATOR) \"" <<
    // v_geodes[position]/*->getParent(0)*/->getName() << "\"" << std::endl;
#else
    _osgView->getCameraManipulator()->setNode(v_geodes[position]);
    // std::cout << "centerMesh \"" << v_geodes[position]->getName() << "\"" <<
    // std::endl;
#endif
    _osgView->getCameraManipulator()->computeHomePosition();
    _osgView->home();


    //FEVV::Debug::print_osg_tree_from_node(v_geodes[position]->getParent(0)->getParent(0));
  }
}


template< typename HalfedgeGraph >
osg::Matrix
FEVV::SimpleViewer< HalfedgeGraph >::getTransformMatrixOsg(unsigned int position)
{
  assert(position < v_geodes.size());
  osg::MatrixTransform *grp_MatrixTransform =
      dynamic_cast< osg::MatrixTransform * >(v_geodes[position]->getParent(0));
  assert(grp_MatrixTransform != nullptr);
  osg::Matrix matrix = grp_MatrixTransform->getMatrix();

  return matrix; // 4x4 homogeneous matrix
}

template< typename HalfedgeGraph >
Eigen::Matrix4d
FEVV::SimpleViewer< HalfedgeGraph >::getTransformMatrixEigen(unsigned int position)
{
  osg::Matrix osg_mat = getTransformMatrixOsg(position);

  // convert OSG transform matrix to Eigen matrix
  // transposition needed!
  Eigen::Matrix4d eigen_mat;
  eigen_mat << osg_mat(0, 0), osg_mat(1, 0), osg_mat(2, 0),    osg_mat(3, 0),
               osg_mat(0, 1), osg_mat(1, 1), osg_mat(2, 1),    osg_mat(3, 1),
               osg_mat(0, 2), osg_mat(1, 2), osg_mat(2, 2),    osg_mat(3, 2),
               osg_mat(0, 3), osg_mat(1, 3), osg_mat(2, 3),    osg_mat(3, 3);
  
  //DBG std::cout << "eigen_mat = \n" << eigen_mat << std::endl;

  return eigen_mat; // 4x4 homogeneous matrix
}

template< typename HalfedgeGraph >
void
FEVV::SimpleViewer< HalfedgeGraph >::resetTransformMatrix(unsigned int position)
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
FEVV::SimpleViewer< HalfedgeGraph >::draw_or_redraw_mesh(
    /*const */ HalfedgeGraph *_g,
    /*const */ PMapsContainer *_pmaps,
    bool _redraw,
    bool _recomputeNT_if_redraw,
    std::string _mesh_filename,
    float _step)
{
  // TODO-elo  fix drawMesh() constness to fix 'mesh' parameter constness

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

  VertexNormalMap *v_nm_ptr = nullptr;
  FaceNormalMap *f_nm_ptr = nullptr;
  VertexColorMap *v_cm_ptr = nullptr;
  FaceColorMap *f_cm_ptr = nullptr;
  VertexUVMap *v_uvm_ptr = nullptr;
  HalfedgeUVMap *h_uvm_ptr = nullptr;

  // textures stuff
  int tex_type = NO_TEXCOORDS;
  std::string texture_filename = "";
  // TODO-elo  remove the 'tex_type' variable ;
  //          the texture type can be deduced of the property maps
  //          that are read from the file ;

  auto pm = get(boost::vertex_point, *_g);

  if(!_redraw)
  {
    m_step = _step;

    // draw mesh
    drawMesh(_g,
             _pmaps,
             &pm,                           /*point map*/
             v_nm_ptr,                      /*vertex-normal map*/
             f_nm_ptr,                      /*face-normal map*/
             v_cm_ptr,                      /*vertex-color map*/
             f_cm_ptr,                      /*face-color map*/
             v_uvm_ptr,                     /*vertex-uv map*/
             h_uvm_ptr,                     /*halfedge-uv map*/
             /* remove */ tex_type,         /*texture type*/
             /* remove */ texture_filename, /*texture filename*/
             _mesh_filename                 /*mesh filename*/
    );

    centerMesh(_g);
  }
  else
  {
    m_recomputeNT_if_redraw = _recomputeNT_if_redraw;

    // redraw mesh
    redrawMesh(true, // here we flush (first param = true) because we add/change
                     // something
               _g,
               _pmaps,
               &pm,                           /*point map*/
               v_nm_ptr,                      /*vertex-normal map*/
               f_nm_ptr,                      /*face-normal map*/
               v_cm_ptr,                      /*vertex-color map*/
               f_cm_ptr,                      /*face-color map*/
               v_uvm_ptr,                     /*vertex-uv map*/
               h_uvm_ptr,                     /*halfedge-uv map*/
               /* remove */ tex_type,         /*texture type*/
               /* remove */ texture_filename, /*texture filename*/
               _mesh_filename                 /*mesh filename*/
    );
  }
}

template< typename HalfedgeGraph >
void
FEVV::SimpleViewer< HalfedgeGraph >::activate_time_mode()
{
  SimpleWindow *sw = static_cast< SimpleWindow * >(
      getWindow()); // here static_cast instead of dynamic_cast only for OSX and
                    // because of plugins... don't understand why...

  sw->activate_time_mode();
}

template< typename HalfedgeGraph >
void
FEVV::SimpleViewer< HalfedgeGraph >::activate_space_mode()
{
  SimpleWindow *sw = static_cast< SimpleWindow * >(
      getWindow()); // here static_cast instead of dynamic_cast only for OSX and
                    // because of plugins... don't understand why...
  
  sw->activate_space_mode();
}

template< typename HalfedgeGraph >
void
FEVV::SimpleViewer< HalfedgeGraph >::setNodeSelected(osg::Node *_geode,
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

    // Assert::check( position != v_meshIsSelected.size(), "can't found this
    // geode", "SimpleViewer::setNodeSelected" ); // @FIXME @jlevallois :
    // disabled since PCL and HEG are in the same file.
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
}

template< typename HalfedgeGraph >
bool
FEVV::SimpleViewer< HalfedgeGraph >::isNodeSelected(osg::Node *_geode)
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

    // Assert::check( position != v_meshIsSelected.size(), "can't found this
    // geode", "SimpleViewer::isNodeSelected" ); // @FIXME @jlevallois :
    // disabled since PCL and HEG are in the same file.
  }

  return false;
}
