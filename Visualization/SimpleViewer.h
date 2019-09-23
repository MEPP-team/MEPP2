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

#include <osgViewer/ViewerEventHandlers>

#include "Base/Color.hpp"
#include "Visualization/BaseViewerOSG.h"

#include "Base/Texture.h"

#include "Visualization/Helpers/OSGHelpers.h"
#include "Visualization/OSG/Visitor/DataVisitor.h"

#include <osg/Material>

#include <osg/Texture2D>
#include <osgDB/ReadFile>

// Generic mesh iterators
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include "FEVV/Wrappings/Graph_traits_extension.h"
#include "FEVV/Wrappings/Geometry_traits.h"
#include <CGAL/boost/graph/properties.h> // included in External folder

#include <iostream>
#include <map>
#include <utility>

#include <Eigen/Dense> // for Eigen::matrix

#ifdef FEVV_USE_CGAL
#include "FEVV/Wrappings/properties_polyhedron_3.h"
#include "FEVV/Wrappings/properties_surface_mesh.h"
#include "FEVV/Wrappings/properties_linear_cell_complex.h"
#include "FEVV/Wrappings/properties_cgal_point_set.h"
#include "FEVV/Wrappings/Graph_traits_extension_cgal_polyhedron_3.h"
#include "FEVV/Wrappings/Graph_traits_extension_cgal_surface_mesh.h"
#include "FEVV/Wrappings/Graph_traits_extension_cgal_linear_cell_complex.h"
#include "FEVV/Wrappings/Graph_traits_extension_cgal_point_set.h"
#endif

#ifdef FEVV_USE_OPENMESH
#include "FEVV/Wrappings/properties_openmesh.h"
#include "FEVV/Wrappings/Graph_traits_extension_openmesh.h"
#endif

#ifdef FEVV_USE_AIF
#include "FEVV/Wrappings/properties_aif.h"
#include "FEVV/Wrappings/Graph_traits_extension_aif.h"
#endif

#ifdef FEVV_USE_PCL
#include "FEVV/Wrappings/properties_pcl_point_cloud.h"
#include "FEVV/Wrappings/Graph_traits_extension_pcl_point_cloud.h"
#endif

namespace FEVV {


/**
 * \brief  Functions to retrieve the name of the datastructure
 *         according to the mesh type.
 */
#ifdef FEVV_USE_CGAL
  inline
  std::string getDatastructureName(FEVV::MeshPolyhedron* m)
  {
    return "POLYHEDRON";
  }

  inline
  std::string getDatastructureName(FEVV::MeshSurface* m)
  {
    return "SURFACEMESH";
  }

  inline
  std::string getDatastructureName(FEVV::MeshLCC* m)
  {
    return "LCC";
  }

  inline
  std::string getDatastructureName(FEVV::CGALPointSet* m)
  {
    return "CGALPOINTSET";
  }
#endif //FEVV_USE_CGAL

#ifdef FEVV_USE_OPENMESH
  inline
  std::string getDatastructureName(FEVV::MeshOpenMesh* m)
  {
    return "OPENMESH";
  }
#endif //FEVV_USE_OPENMESH

#ifdef FEVV_USE_AIF
  inline
  std::string getDatastructureName(FEVV::MeshAIF* m)
  {
    return "AIF";
  }
#endif //FEVV_USE_AIF

#ifdef FEVV_USE_PCL
  inline
  std::string getDatastructureName(FEVV::PCLPointCloud* m)
  {
    return "PCLPOINTCLOUD";
  }
#endif //FEVV_USE_PCL


/**
 * \brief  A container to store pointers over meshes of mixed types.
 */
class MixedMeshesVector
{
public:
  typedef  std::pair<void*, std::string>  VoidMeshPtr;

  void push_back(VoidMeshPtr &ptr)
  {
    storage.push_back(ptr);
  }

  template< typename MeshT >
  void push_back(MeshT *m)
  {
    auto pair = std::make_pair(static_cast< void* >(m),
                               getDatastructureName(m));
    storage.push_back(pair);
    std::cout << "[SimpleViewer] mesh " << pair.first
              << " stored with datastructure " << pair.second
              << std::endl;
  }

  std::size_t size(void)
  {
    return storage.size();
  }

  VoidMeshPtr&  operator[](std::size_t pos)
  {
    return storage[pos];
  }

protected:
  std::vector< VoidMeshPtr > storage;
};


/**
 * \class SimpleViewer
 * \brief SimpleViewer is a specialization of osgViewer::CompositeViewer.
 * This class is a widget where we are able to draw objects using
 * OpenSceneGraph.
 */
class SimpleViewer : public BaseViewerOSG
{
public:
  /**
   * Constructor.
   */
  SimpleViewer();

  ~SimpleViewer();

  void init() override;

  bool isInit() const override;

  bool isValid() const override;

  /**
   *
   * @note Actually, it returns true only if at least one mesh is selected.
   */
  // bool isSelected() const override; MT -> JL : why this function ?

  // void attach( Widget* _widget );

  /**
   * Change background color of the scene.
   *
   * @param[in] _color the color to use.
   **/
  bool changeBackgroundColor(const Color &_color) override;

  /**
   * Export the current view of the scene into a screenshot (JPEG).
   *
   * @note Can not be const due to OSG.
   *
   * @param[in] _name name of the file to export (without extension).
   **/
  bool saveScreenshot(const std::string &_name) override;

  /**
   * Add a geode to the scene.
   *
   * @note A geode is a Geometry Node.
   *
   * @param[in] _geode Pointer to a geode.
   **/
  void addModel(Model *_geode) override;

  /**
   * Add a group to the scene.
   *
   * @note A group is a set of geodes.
   *
   * @param[in] _group Pointer to a group of geode.
   **/
  void addGroup(Group *_group) override;

  void setNodeSelected(osg::Node *_geode, bool isSelected) override;
  bool isNodeSelected(osg::Node *_geode) override;

  /**
   * Returns all selected geodes.
   *
   * @note If none geode is selected, the return will be empty.
   *
   * @return a std::vector of pointers of all selected geodes.
   */
  std::vector< osg::Geode * > getSelectedGeodes();

  /**
   * Returns all geodes.
   *
   * @note If there is no mesh, the return will be empty.
   *
   * @return a std::vector of pointers of all geodes.
   */
  std::vector< osg::Geode * > getGeodes();

  /**
   * Returns all selected meshes.
   *
   * @note If none mesh is selected, the return will be empty.
   *
   * @return a std::vector of pointers of all selected meshes.
   */
  MixedMeshesVector getSelectedMeshes();

  /**
   * Returns all meshes.
   *
   * @note If there is no mesh, the return will be empty.
   *
   * @return a std::vector of pointers of all meshes.
   */
  MixedMeshesVector getMeshes();

  /**
   * Returns the id of the mesh in v_mixed_meshes array.
   * Returns -1 if the mesh is not in v_mixed_meshes array.
   *
   * @return  index of the mesh in v_mixed_meshes array.
   */
  size_t getMeshId(const void *mesh_ptr);

  /**
   * Returns all selected meshes names.
   *
   * @note If none mesh is selected, the return will be empty.
   *
   * @return a std::vector of all selected meshes names.
   */
  std::vector< std::string > getSelectedMeshesNames();

  /**
   * Returns all meshes names.
   *
   * @note If there is no mesh, the return will be empty.
   *
   * @return a std::vector of all meshes names.
   */
  std::vector< std::string > getMeshesNames();

  /**
   * Returns all selected properties maps.
   *
   * @note If none mesh is selected, the return will be empty.
   *
   * @return a std::vector of pointers of all selected properties maps.
   */
  std::vector< PMapsContainer * > getSelected_properties_maps();

  /**
   * Returns all properties maps.
   *
   * @note If there is no properties map, the return will be empty.
   *
   * @return a std::vector of pointers of all properties maps.
   */
  std::vector< PMapsContainer * > get_properties_maps();


  std::vector< osg::Group * > getSelectedDraggers1();
  std::vector< osg::Group * > getDraggers1();

  std::vector< osg::Group * > getSelectedDraggers2();
  std::vector< osg::Group * > getDraggers2();


#if 0 //TODO-elo-rm-?-ask_MTO
  /**
   * Returns the mesh at a given position.
   *
   * @note If there is no mesh at the given positon, an assert will be raised.
   *
   * @param[in] _position A given position (lower than v_meshes.size()).
   *
   * @return a pointer of the mesh at a given position.
   */
  HalfedgeGraph *getMesh(unsigned int _position);
#endif

  DataModelVector *getDataModel() override;


  /**
   * Returns the transformation matrix of the mesh at a given position.
   *
   * @note If there is no mesh at the given positon, an assert will be raised.
   *
   * @param[in] position A given position (lower than v_meshes.size()).
   *
   * @return a 4x4 homogeneous matrix.
   */
  osg::Matrix getTransformMatrixOsg(unsigned int position);
  Eigen::Matrix4d getTransformMatrixEigen(unsigned int position);

  /**
   * Reset the transformation matrix of the mesh at a given position.
   *
   * @note If there is no mesh at the given positon, an assert will be raised.
   *
   * @param[in] position A given position (lower than v_meshes.size()).
   */
  void resetTransformMatrix(unsigned int position);


  /**
   * Draw mesh into the scene.
   *
   * @note The mesh must be a model of HalfedgeGraph concept.
   *
   * @tparam      PointMap    a class of points map.
   *                          Associate a coordinate to an
   *                          halfedge_descriptor.
   *
   * @param[in]   _g          a mesh (model of HalfedgeGraph concept).
   * @param[in]   _pm         a point map.
   **/
  template< typename HalfedgeGraph, typename PointMap >
  osg::ref_ptr< osg::Group >
  createMesh(HalfedgeGraph *_g,
             PMapsContainer *_pmaps,
             PointMap *_pm,
             std::string _mesh_file = std::string(""),
             osg::ref_ptr< osg::Group > _group = new osg::Group);

  template< typename HalfedgeGraph, typename PointMap >
  void drawMesh(HalfedgeGraph *_g,
                PMapsContainer *_pmaps,
                PointMap *_pm,
                std::string _mesh_file = std::string(""));

  template< typename HalfedgeGraph, typename PointMap >
  void redrawMesh(HalfedgeGraph *_g,
                  PMapsContainer *_pmaps,
                  PointMap *_pm = nullptr,
                  std::string _mesh_file = std::string(""));

  template< typename HalfedgeGraph >
  void centerMesh(HalfedgeGraph *_g);

  template< typename HalfedgeGraph >
  void draw_or_redraw_mesh(/*const */ HalfedgeGraph *_g,
                           /*const */ PMapsContainer *_pmaps,
                           bool _redraw = false,
                           bool _recomputeNT_if_redraw = false,
                           std::string _mesh_filename = std::string(""),
                           bool _recreateOSGobj_if_redraw = true,
                           float _step = 0.);

  void activate_time_mode();
  void activate_space_mode();

  void updateSWModelList();

  osg::Node *addDraggersToScene(osg::Node *scene,
                                const std::string &nameDrag1,
                                float fScaleDrag1,
                                char keyDrag1,
                                const std::string &nameDrag2,
                                float fScaleDrag2,
                                char keyDrag2);

protected:
  /**
   * Draw mesh into the scene.
   *
   * @note The mesh must be a model of HalfedgeGraph concept.
   *
   * @tparam      PointMap    a class of points map.
   *                          Associate a coordinate to an
   *                          halfedge_descriptor.
   *
   * @param[in]   _g          a mesh (model of HalfedgeGraph concept).
   * @param[in]   _pm         a point map.
   **/
  template< typename HalfedgeGraph, typename PointMap >
  void internal_createMesh(
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
      std::string _mesh_file = std::string(""));
  /**
   * Draw mesh into the scene.
   *
   * @note The mesh must be a model of HalfedgeGraph concept.
   *
   * @tparam      PointMap    a class of points map.
   *                          Associate a coordinate to an
   *                          halfedge_descriptor.
   *
   * @param[in]   _g          a mesh (model of HalfedgeGraph concept).
   * @param[in]   _pm         a point map.
   **/
  template< typename HalfedgeGraph, typename PointMap >
  osg::Geode *internal_createMesh(
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
      std::string _mesh_file = std::string(""));

  /**
   * Draw point cloud into the scene.
   *
   * @note The point cloud must be a model of PointCloud concept.
   *
   * @tparam      PointMap    a class of points map.
   *                          Associate a coordinate to an
   *                          vertex_descriptor.
   *
   * @param[in]   _g          a point cloud (model of PointCloud concept).
   * @param[in]   _pm         a point map.
   **/
  template< typename PointCloud, typename PointMap >
  void internal_createMesh_pointcloud(
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
      std::string _mesh_file = std::string(""));

private:
  /**
   * Loads a mesh into the scene, using shaders.
   * Assigns a shader to each mesh, according to its material type (standard -
   *Blinn-Phong/PBR - Cook-Torrance).
   *
   * @note The mesh must be a model of HalfedgeGraph concept.
   *
   * @tparam      VertexNormalMap    a class of vertex normals map.
   *                                 Associates a normal to a
   *                                 vertex_descriptor.
   * @tparam      VertexTangentMap    a class of vertex tangents map.
   *                                  Associates a tangent to an
   *                                  vertex_descriptor.
   * @tparam      VertexColorMap    a class of vertex colors map.
   *                                Associates a coordinate to a
   *                                vertex_descriptor.
   * @tparam      FaceMaterialMap    a class of face materials map.
   *                                 Associates a material to an
   *                                 material index.
   *
   * @param[in]   _geode                 the geode which we want to load.
   * @param[in]   _g                     the halfedge graph.
   * @param[in]   _geometries            the geometries which we load into the
   *geode.
   * @param[in]   _geometries_edges      the superimposed edges geometries which
   *we load into the geode.
   * @param[in]   _geometries_vertices   the superimposed vertices geometries
   *which we load into the geode.
   * @param[in]   _vertexArrays          the vertices positions we load into the
   *geode.
   * @param[in]   _vertexArrays_edges    the superimposed edges positions we
   *load into the geode.
   * @param[in]   _vertexArrays_vertices the superimposed vertices positions we
   *load into the geode.
   * @param[in]   _normalsArrays         the normals we load into the geode.
   * @param[in]   _tangentsArrays        the tangents we load into the geode.
   * @param[in]   _texcoordsArrays       the texture coordinates we load into
   *the geode.
   * @param[in]   _colorsArrays          the colors we load into the geode.
   * @param[in]   _colorsArrays_edges    the superimposed edges colors we load
   *into the geode.
   * @param[in]   _colorsArrays_vertices the superimposed vertices colors we
   *load into the geode.
   * @param[in]   _m_mm_size             the number of materials.
   * @param[in]   _texture_type          the texture coordinates setup (by
   *vertex, by halfedge, or none).
   * @param[in]   _vt_nm                 the vertex normals property map.
   * @param[in]   _vt_tm                 the vertex tangents property map.
   * @param[in]   _vt_cm                 the vertex colors property map.
   * @param[in]   _m_mm                  the face materials property map.
   **/
  template< typename HalfedgeGraph,
            typename VertexNormalMap,
            typename VertexTangentMap,
            typename VertexColorMap,
            typename FaceColorMap,
            typename VertexUVMap,
            typename HalfedgeUVMap,
            typename FaceMaterialMap >
  void internal_loadShadedMesh(
      osg::Geode *_geode,
      HalfedgeGraph *_g,
      const std::vector< osg::ref_ptr< osg::Geometry > > &_geometries,
      const std::vector< osg::ref_ptr< osg::Geometry > > &_geometries_edges,
      const std::vector< osg::ref_ptr< osg::Geometry > > &_geometries_vertices,
      const std::vector< osg::ref_ptr< osg::Geometry > > &_geometries_normals,
      const std::vector< osg::ref_ptr< osg::Vec3Array > > &_vertexArrays,
      const std::vector< osg::ref_ptr< osg::Vec3Array > > &_vertexArrays_edges,
      const std::vector< osg::ref_ptr< osg::Vec3Array > > &_vertexArrays_vertices,
      const std::vector< osg::ref_ptr< osg::Vec3Array > > &_vertexArrays_normals,
      const std::vector< osg::ref_ptr< osg::Vec3Array > > &_normalsArrays,
      const std::vector< osg::ref_ptr< osg::Vec3Array > > &_normalsArrays_edges,
      const std::vector< osg::ref_ptr< osg::Vec3Array > > &_normalsArrays_vertices,
      const std::vector< osg::ref_ptr< osg::Vec3Array > > &_tangentsArrays,
      const std::vector< osg::ref_ptr< osg::Vec2Array > > &_texcoordsArrays,
      const std::vector< osg::ref_ptr< osg::Vec4Array > > &_colorsArrays,
      const std::vector< osg::ref_ptr< osg::Vec4Array > > &_colorsArrays_edges,
      const std::vector< osg::ref_ptr< osg::Vec4Array > > &_colorsArrays_vertices,
      const std::vector< osg::ref_ptr< osg::Vec4Array > > &_colorsArrays_normals,
      std::size_t _m_mm_size,
      VertexNormalMap *_vt_nm,
      VertexTangentMap *_vt_tm,
      VertexColorMap *_vt_cm,
      FaceColorMap *_f_cm,
      VertexUVMap *_vt_uv_m,
      HalfedgeUVMap *_het_uv_m,
      FaceMaterialMap *_m_mm);

  /**
   * Loads a mesh into the scene using old rendering mode.
   *
   * @note The mesh must be a model of HalfedgeGraph concept.
   *
   * @tparam      VertexNormalMap    a class of vertex normals map.
   *                                 Associates a normal to a
   *                                 vertex_descriptor.
   * @tparam      VertexColorMap    a class of vertex colors map.
   *                                Associates a coordinate to a
   *                                vertex_descriptor.
   * @tparam      FaceColorMap    a class of face colors map.
   *                              Associates a coordinate to a
   *                              material index.
   * @tparam      VertexUVMap    a class of vertex texture coordinates map.
   *                             Associates texture coordinates to a
   *                             vertex_descriptor.
   * @tparam      HalfedgeUVMap    a class of halfedge texture coordinates map.
   *                               Associates texture coordinates to an
   *                               halfedge_descriptor.
   * @tparam      FaceMaterialMap    a class of face materials map.
   *                                 Associates a material to an
   *                                 material index.
   *
   * @param[in]   _geode                 the geode which we want to load.
   * @param[in]   _g                     the halfedge graph.
   * @param[in]   _geometries            the geometries which we load into the
   *geode.
   * @param[in]   _geometries_edges      the superimposed edges geometries which
   *we load into the geode.
   * @param[in]   _geometries_vertices   the superimposed vertices geometries
   *which we load into the geode.
   * @param[in]   _vertexArrays          the vertices positions we load into the
   *geode.
   * @param[in]   _vertexArrays_edges    the superimposed edges positions we
   *load into the geode.
   * @param[in]   _vertexArrays_vertices the superimposed vertices positions we
   *load into the geode.
   * @param[in]   _normalsArrays         the normals we load into the geode.
   * @param[in]   _texcoordsArrays       the texture coordinates we load into
   *the geode.
   * @param[in]   _colorsArrays          the colors we load into the geode.
   * @param[in]   _colorsArrays_edges    the superimposed edges colors we load
   *into the geode.
   * @param[in]   _colorsArrays_vertices the superimposed vertices colors we
   *load into the geode.
   * @param[in]   _m_mm_size             the number of materials.
   * @param[in]   _texture_type          the texture coordinates setup (by
   *vertex, by halfedge, or none).
   * @param[in]   _vt_nm                 the vertex normals property map.
   * @param[in]   _vt_cm                 the vertex colors property map.
   * @param[in]   _f_cm                  the face colors property map.
   * @param[in]   _vt_uv_m               the vertex texture coordinates property
   *map.
   * @param[in]   _het_uv_m              the halfedge texture coordinates
   *property map.
   * @param[in]   _m_mm                  the face materials property map.
   **/
  template< typename HalfedgeGraph,
            typename VertexNormalMap,
            typename VertexColorMap,
            typename FaceColorMap,
            typename VertexUVMap,
            typename HalfedgeUVMap,
            typename FaceMaterialMap >
  void internal_loadLegacyMesh(
      osg::Geode *_geode,
      HalfedgeGraph *_g,
      const std::vector< osg::ref_ptr< osg::Geometry > > &_geometries,
      const std::vector< osg::ref_ptr< osg::Geometry > > &_geometries_edges,
      const std::vector< osg::ref_ptr< osg::Geometry > > &_geometries_vertices,
      const std::vector< osg::ref_ptr< osg::Geometry > > &_geometries_normals,
      const std::vector< osg::ref_ptr< osg::Vec3Array > > &_vertexArrays,
      const std::vector< osg::ref_ptr< osg::Vec3Array > > &_vertexArrays_edges,
      const std::vector< osg::ref_ptr< osg::Vec3Array > > &_vertexArrays_vertices,
      const std::vector< osg::ref_ptr< osg::Vec3Array > > &_vertexArrays_normals,
      const std::vector< osg::ref_ptr< osg::Vec3Array > > &_normalsArrays,
      const std::vector< osg::ref_ptr< osg::Vec3Array > > &_normalsArrays_edges,
      const std::vector< osg::ref_ptr< osg::Vec3Array > > &_normalsArrays_vertices,
      const std::vector< osg::ref_ptr< osg::Vec2Array > > &_texcoordsArrays,
      const std::vector< osg::ref_ptr< osg::Vec4Array > > &_colorsArrays,
      const std::vector< osg::ref_ptr< osg::Vec4Array > > &_colorsArrays_edges,
      const std::vector< osg::ref_ptr< osg::Vec4Array > > &_colorsArrays_vertices,
      const std::vector< osg::ref_ptr< osg::Vec4Array > > &_colorsArrays_normals,
      std::size_t _m_mm_size,
      int _texture_type,
      VertexNormalMap *_vt_nm,
      VertexColorMap *_vt_cm,
      FaceColorMap *_f_cm,
      VertexUVMap *_vt_uv_m,
      HalfedgeUVMap *_het_uv_m,
      FaceMaterialMap *_m_mm);

  // List of lights used in the scene
  std::vector< osg::ref_ptr< osg::Light > > lights;

protected:
  MixedMeshesVector v_mixed_meshes;

  std::vector< std::string > v_meshes_names;
  std::vector< PMapsContainer * > v_properties_maps;
  std::vector< osg::Group * > v_draggers1;
  std::vector< osg::Group * > v_draggers2;
  std::vector< bool > v_meshIsSelected;
  std::vector< osg::Geode * > v_geodes;

  std::vector< std::vector< osg::ref_ptr< osg::Geometry > > > v_geometries,
      v_geometriesL, v_geometriesP, v_geometries_edges, v_geometries_vertices, v_geometries_normals;
  std::vector< std::vector< osg::ref_ptr< osg::Vec3Array > > > v_vertexArrays,
      v_vertexArrays_edges, v_vertexArrays_vertices, v_vertexArrays_normals;
  std::vector< std::vector< osg::ref_ptr< osg::Vec3Array > > > v_normalsArrays,
      v_normalsArraysF, v_normalsArrays_edges, v_normalsArrays_vertices;
  std::vector< std::vector< osg::ref_ptr< osg::Vec3Array > > > v_tangentsArrays;
  std::vector< std::vector< osg::ref_ptr< osg::Vec4Array > > > v_colorsArrays,
      v_colorsArrays_edges, v_colorsArrays_vertices, v_colorsArrays_normals;
  std::vector< std::vector< osg::ref_ptr< osg::Vec2Array > > >
      v_texcoordsArrays;

public:
  osg::ref_ptr< osg::Group > gizmo;
  osg::ref_ptr< osg::Group > grid;

  osg::ref_ptr< osg::Group > hud;
  osg::ref_ptr< osgText::Text > hudText;

  int i_time;
  int current_i_time;
};

} // namespace FEVV

#include "Visualization/MeshLoading.inl"
#include "Visualization/SimpleViewer.inl"

#ifdef FEVV_USE_CGAL
#include "Visualization/SimpleViewerCGALPointSet.inl"
#endif

#ifdef FEVV_USE_PCL
#include "Visualization/SimpleViewerPCLPointCloud.inl"
#endif