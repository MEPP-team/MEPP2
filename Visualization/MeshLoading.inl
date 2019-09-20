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

#include <osg/ShadeModel>
#include <osgUtil/Optimizer>

#include <memory>

#include "Visualization/SimpleWindow.h"

// TODO : REVOIR les osg::StateAttribute::OVERRIDE et les vérifier !!!

struct CameraPosCallback : public osg::Uniform::Callback
{
  explicit CameraPosCallback(osg::Camera *camera) : camera{camera} {}

  virtual void operator()(osg::Uniform *uniform,
                          osg::NodeVisitor *nodeVisitor) override
  {
    uniform->set(osg::Vec3(camera->getInverseViewMatrix().getTrans()));
  }

  osg::Camera *camera;
};

struct CameraTargetCallback : public osg::Uniform::Callback
{
  explicit CameraTargetCallback(osg::Geode *geode) : geode{geode} {}

  virtual void operator()(osg::Uniform *uniform,
                          osg::NodeVisitor *nodeVisitor) override
  {
    uniform->set(geode->getBound().center());

    // The bounding box being never recomputed, this line should work although
    // less accurate
    // uniform->set(osg::Vec3(osg::computeLocalToWorld(nodeVisitor->getNodePath()).getTrans()));
  }

  osg::Geode *geode;
};

struct ModelMatrixCallback : public osg::Uniform::Callback
{
  virtual void operator()(osg::Uniform *uniform,
                          osg::NodeVisitor *nodeVisitor) override
  {
    uniform->set(osg::computeLocalToWorld(nodeVisitor->getNodePath()));
  }
};

struct ViewMatrixCallback : public osg::Uniform::Callback
{
  explicit ViewMatrixCallback(osg::Camera *camera) : camera{camera} {}

  virtual void operator()(osg::Uniform *uniform,
                          osg::NodeVisitor *nodeVisitor) override
  {
    uniform->set(osg::Matrix(camera->getViewMatrix()));
  }

  osg::Camera *camera;
};

struct InverseViewMatrixCallback : public osg::Uniform::Callback
{
  explicit InverseViewMatrixCallback(osg::Camera *camera) : camera{camera} {}

  virtual void operator()(osg::Uniform *uniform,
                          osg::NodeVisitor *nodeVisitor) override
  {
    uniform->set(osg::Matrix(camera->getInverseViewMatrix()));
  }

  osg::Camera *camera;
};

struct ProjectionMatrixCallback : public osg::Uniform::Callback
{
  explicit ProjectionMatrixCallback(osg::Camera *camera) : camera{camera} {}

  virtual void operator()(osg::Uniform *uniform,
                          osg::NodeVisitor *nodeVisitor) override
  {
    uniform->set(osg::Matrix(camera->getProjectionMatrix()));
  }

  osg::Camera *camera;
};

struct MVPMatrixCallback : public osg::Uniform::Callback
{
  explicit MVPMatrixCallback(osg::Camera *camera) : camera{camera} {}

  virtual void operator()(osg::Uniform *uniform,
                          osg::NodeVisitor *nodeVisitor) override
  {
    const auto viewMatrix = camera->getViewMatrix();
    const auto modelMatrix =
        osg::computeLocalToWorld(nodeVisitor->getNodePath());
    const auto mvpMatrix =
        modelMatrix * viewMatrix * camera->getProjectionMatrix();

    uniform->set(mvpMatrix);
  }

  osg::Camera *camera;
};

inline std::string
recoverBaseDirectory()
{
  auto baseDir = QDir(QApplication::applicationDirPath());

#if defined(Q_OS_MAC)
  if(baseDir.dirName() == "MacOS")
  {
    baseDir.cdUp();
    baseDir.cdUp();
    baseDir.cdUp();
  }
#endif

  return baseDir.absolutePath().toStdString();
}

inline osg::ref_ptr< osg::Texture2D >
createDefaultTexture(uint8_t value = 255) // white by default
{
  osg::ref_ptr< osg::Image > img = new osg::Image();
  img->allocateImage(1, 1, 1, GL_RGB, GL_UNSIGNED_BYTE);

  const osg::Vec3ub color(value, value, value);
  const auto imgData = reinterpret_cast< osg::Vec3ub * >(img->data());
  *imgData = color;

  osg::ref_ptr< osg::Texture2D > texture = new osg::Texture2D();
  texture->setWrap(osg::Texture2D::WRAP_S, osg::Texture2D::CLAMP_TO_EDGE);
  texture->setWrap(osg::Texture2D::WRAP_T, osg::Texture2D::CLAMP_TO_EDGE);
  texture->setFilter(osg::Texture2D::MIN_FILTER, osg::Texture2D::NEAREST);
  texture->setFilter(osg::Texture2D::MAG_FILTER, osg::Texture2D::NEAREST);
  texture->setImage(img);

  return texture;
}

inline osg::ref_ptr< osg::Texture2D >
createTexture(const std::string &filename)
{
  osg::ref_ptr< osg::Texture2D > texture;

  // load an image by reading a file
  const osg::ref_ptr< osg::Image > image = osgDB::readImageFile(filename);
  if(image)
  {
    // std::cout << "[MeshLoading] Texture file '" << _texture_file << "'
    // loaded." << std::endl;

    texture = new osg::Texture2D;

    // protect from being optimized away as static state
    texture->setDataVariance(osg::Object::STATIC); // previously DYNAMIC

    // assign the texture to the image we read from file
    texture->setImage(image);

    texture->setWrap(osg::Texture2D::WRAP_S, osg::Texture2D::REPEAT);
    texture->setWrap(osg::Texture2D::WRAP_T, osg::Texture2D::REPEAT);
    texture->setFilter(osg::Texture2D::MIN_FILTER,
                       osg::Texture2D::LINEAR_MIPMAP_LINEAR);
    texture->setFilter(osg::Texture2D::MAG_FILTER, osg::Texture2D::LINEAR);
  }
  else
  {
    std::cerr << "-> [MeshLoading] Couldn't find texture file, creating "
                 "default texture."
              << std::endl;
    texture = createDefaultTexture();
  }

  return texture;
}

inline void
loadMaterialStandard(const osg::ref_ptr< osg::Geometry > &geometry,
                     const FEVV::Types::Material &material)
{
  uint8_t map_index = 0;

  geometry->getOrCreateStateSet()->addUniform(
      new osg::Uniform("uniMaterial.ambientMap", 0));
  geometry->getOrCreateStateSet()->addUniform(
      new osg::Uniform("uniMaterial.diffuseMap", 1));
  geometry->getOrCreateStateSet()->addUniform(
      new osg::Uniform("uniMaterial.specularMap", 2));
  geometry->getOrCreateStateSet()->addUniform(
      new osg::Uniform("uniMaterial.emissiveMap", 3));
  geometry->getOrCreateStateSet()->addUniform(
      new osg::Uniform("uniMaterial.transparencyMap", 4));
  geometry->getOrCreateStateSet()->addUniform(
      new osg::Uniform("uniMaterial.bumpMap", 5));

  for(const auto &texture_filename : {material.ambient_texture_filename,
                                      material.diffuse_texture_filename,
                                      material.specular_texture_filename,
                                      material.emissive_texture_filename,
                                      material.transparency_texture_filename,
                                      material.normal_map_filename})
  {
    osg::ref_ptr< osg::Texture2D > texture;

    if(!texture_filename.empty())
      texture = createTexture(texture_filename);
    else
      texture = createDefaultTexture(
          (map_index == 3 ? 0
                          : 255)); // Little hack here: for emissive map, create
                                   // a black texture instead of a white

    geometry->getOrCreateStateSet()->setTextureAttribute(
        map_index,
        texture,
        osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
    ++map_index;
  }
}

inline void
loadMaterialPBR(const osg::ref_ptr< osg::Geometry > &geometry,
                const FEVV::Types::Material &material)
{
  uint8_t map_index = 0;

  geometry->getOrCreateStateSet()->addUniform(
      new osg::Uniform("uniMaterial.albedoMap", 0));
  geometry->getOrCreateStateSet()->addUniform(
      new osg::Uniform("uniMaterial.normalMap", 1));
  geometry->getOrCreateStateSet()->addUniform(
      new osg::Uniform("uniMaterial.metallicMap", 2));
  geometry->getOrCreateStateSet()->addUniform(
      new osg::Uniform("uniMaterial.roughnessMap", 3));
  geometry->getOrCreateStateSet()->addUniform(
      new osg::Uniform("uniMaterial.emissiveMap", 4));
  geometry->getOrCreateStateSet()->addUniform(
      new osg::Uniform("uniMaterial.ambientOcclusionMap", 5));

  for(const auto &texture_filename : {material.diffuse_texture_filename,
                                      material.normal_map_filename,
                                      material.metallic_map_filename,
                                      material.roughness_map_filename,
                                      material.emissive_texture_filename,
                                      material.ambient_texture_filename})
  {
    osg::ref_ptr< osg::Texture2D > texture;

    if(!texture_filename.empty())
      texture = createTexture(texture_filename);
    else
      texture = createDefaultTexture(
          (map_index == 4 ? 0
                          : 255)); // Little hack here: for emissive map, create
                                   // a black texture instead of a white

    geometry->getOrCreateStateSet()->setTextureAttribute(
        map_index,
        texture,
        osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
    ++map_index;
  }
}

template< typename VertexNormalMap,
          typename VertexTangentMap,
          typename FaceMaterialMap >
inline osg::ref_ptr< osg::Program >
loadTexturedMesh(
    const std::vector< osg::ref_ptr< osg::Geometry > > &_geometries,
    const std::vector< osg::ref_ptr< osg::Vec3Array > > &_vertexArrays,
    const std::vector< osg::ref_ptr< osg::Vec3Array > > &_normalsArrays,
    const std::vector< osg::ref_ptr< osg::Vec3Array > > &_tangentsArrays,
    const std::vector< osg::ref_ptr< osg::Vec2Array > > &_texcoordsArrays,
    bool _useSmoothShading,
    VertexNormalMap *_vt_nm,
    VertexTangentMap *_vt_tm,
    FaceMaterialMap *_m_mm,
    std::size_t unit_ii)
{
  const auto &shadersDirectory = recoverBaseDirectory() + "/Shaders/";
  osg::ref_ptr< osg::Program > program = new osg::Program;

  osg::ref_ptr< osg::Shader > vertShader = new osg::Shader(osg::Shader::VERTEX);
  if(!vertShader->loadShaderSourceFromFile(shadersDirectory + "vert.glsl"))
    std::cerr << "[MeshLoading] Could not read VERTEX shader from file"
              << std::endl;
  program->addShader(vertShader);

  // RM: disabled these calls, fixing loading time & FPS problems
  //_geometries[unit_ii]->setUseVertexBufferObjects(true);

  // Sending texcoords to shader
  _geometries[unit_ii]->setVertexAttribArray(
      1, _texcoordsArrays[unit_ii], osg::Array::BIND_PER_VERTEX);

  // Sending tangents to shader (no need if we're rendering in flat shading)
  if(_vt_tm != nullptr && _useSmoothShading)
    _geometries[unit_ii]->setVertexAttribArray(
        3, _tangentsArrays[unit_ii], osg::Array::BIND_PER_VERTEX);

  std::string fragShaderLocation =
      shadersDirectory + "blinn-phong.glsl"; // Shader by default

  const auto &material = get(*_m_mm, unit_ii);
  if(!material.name.empty()) // Normally, 'empty' is not possible here (would be
                             // a colored mesh)
  {
    if(material.type == FEVV::Types::MaterialType::MATERIAL_TYPE_PBR)
    {
      fragShaderLocation = shadersDirectory + "cook-torrance.glsl";

      _geometries[unit_ii]->getOrCreateStateSet()->addUniform(
          new osg::Uniform("uniMaterial.baseColor",
                           osg::Vec3f(material.diffuse_red_component,
                                      material.diffuse_green_component,
                                      material.diffuse_blue_component)));
      _geometries[unit_ii]->getOrCreateStateSet()->addUniform(
          new osg::Uniform("uniMaterial.metallicFactor",
                           static_cast< float >(material.metallic_factor)));
      _geometries[unit_ii]->getOrCreateStateSet()->addUniform(
          new osg::Uniform("uniMaterial.roughnessFactor",
                           static_cast< float >(material.roughness_factor)));
      _geometries[unit_ii]->getOrCreateStateSet()->addUniform(
          new osg::Uniform("uniMaterial.emissive",
                           osg::Vec3f(material.emissive_red_component,
                                      material.emissive_green_component,
                                      material.emissive_blue_component)));

      loadMaterialPBR(_geometries[unit_ii], material);
    }
    else
    {
      _geometries[unit_ii]->getOrCreateStateSet()->addUniform(
          new osg::Uniform("uniMaterial.ambient",
                           osg::Vec3f(material.ambient_red_component,
                                      material.ambient_green_component,
                                      material.ambient_blue_component)));
      _geometries[unit_ii]->getOrCreateStateSet()->addUniform(
          new osg::Uniform("uniMaterial.diffuse",
                           osg::Vec3f(material.diffuse_red_component,
                                      material.diffuse_green_component,
                                      material.diffuse_blue_component)));
      _geometries[unit_ii]->getOrCreateStateSet()->addUniform(
          new osg::Uniform("uniMaterial.specular",
                           osg::Vec3f(material.specular_red_component,
                                      material.specular_green_component,
                                      material.specular_blue_component)));
      _geometries[unit_ii]->getOrCreateStateSet()->addUniform(
          new osg::Uniform("uniMaterial.emissive",
                           osg::Vec3f(material.emissive_red_component,
                                      material.emissive_green_component,
                                      material.emissive_blue_component)));
      _geometries[unit_ii]->getOrCreateStateSet()->addUniform(
          new osg::Uniform("uniMaterial.transparency",
                           static_cast< float >(material.transparency)));

      loadMaterialStandard(_geometries[unit_ii], material);
    }

    if(material.has_normal_map && _useSmoothShading)
      _geometries[unit_ii]->getOrCreateStateSet()->addUniform(
          new osg::Uniform("uniUseNormalMapping", true));
  }

  osg::ref_ptr< osg::Shader > fragShader =
      new osg::Shader(osg::Shader::FRAGMENT);
  if(!fragShader->loadShaderSourceFromFile(fragShaderLocation))
    std::cerr << "[MeshLoading] Could not read FRAGMENT shader from file"
              << std::endl;

  program->addShader(fragShader);

  return program;
}

template< typename VertexNormalMap, typename VertexColorMap >
inline osg::ref_ptr< osg::Program >
loadColoredMesh(
    const std::vector< osg::ref_ptr< osg::Geometry > > &_geometries,
    const std::vector< osg::ref_ptr< osg::Vec3Array > > &_vertexArrays,
    const std::vector< osg::ref_ptr< osg::Vec3Array > > &_normalsArrays,
    const std::vector< osg::ref_ptr< osg::Vec4Array > > &_colorsArrays,
    VertexNormalMap *_vt_nm,
    VertexColorMap *_vt_cm,
    std::size_t unit_ii)
{
  const auto &shadersDirectory = recoverBaseDirectory() + "/Shaders/";
  osg::ref_ptr< osg::Program > program = new osg::Program;

  osg::ref_ptr< osg::Shader > vertShader = new osg::Shader(osg::Shader::VERTEX);
  if(!vertShader->loadShaderSourceFromFile(shadersDirectory + "colorVert.glsl"))
    std::cerr << "[MeshLoading] Could not read VERTEX shader from file"
              << std::endl;
  program->addShader(vertShader);

  osg::ref_ptr< osg::Shader > fragShader =
      new osg::Shader(osg::Shader::FRAGMENT);
  if(!fragShader->loadShaderSourceFromFile(shadersDirectory + "colorFrag.glsl"))
    std::cerr << "[MeshLoading] Could not read FRAGMENT shader from file"
              << std::endl;
  program->addShader(fragShader);

  // RM: disabled these calls, fixing loading time & FPS problems
  //_geometries[unit_ii]->setUseVertexBufferObjects(true);

  // Sending colors to shader
  _geometries[unit_ii]->setVertexAttribArray(
      1,
      _colorsArrays[unit_ii],
      (_vt_cm != nullptr ? osg::Array::BIND_PER_VERTEX
                         : osg::Array::BIND_PER_PRIMITIVE_SET));

  return program;
}

template< typename HalfedgeGraph,
          typename VertexNormalMap,
          typename VertexTangentMap,
          typename VertexColorMap,
          typename FaceColorMap,
          typename VertexUVMap,
          typename HalfedgeUVMap,
          typename FaceMaterialMap >
void
FEVV::SimpleViewer::internal_loadShadedMesh(
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
    FaceMaterialMap *_m_mm)
{
  std::cout << "[MeshLoading] Loading using shaders." << std::endl;

  const auto &shadersDirectory = recoverBaseDirectory() + "/Shaders/";

  std::vector< osg::ref_ptr< osg::Vec4Array > > colorsArrays2;

  size_t unit_ii = 0;
  do
  {
    _geometries[unit_ii]->setStateSet(NULL);   // NEW
    _geometries[unit_ii]->setColorArray(NULL); // NEW

    osg::ref_ptr< osg::Program > program;

    // Sending positions to shader
    _geometries[unit_ii]->setVertexAttribArray(
        0, _vertexArrays[unit_ii], osg::Array::BIND_PER_VERTEX);

    // Sending normals to shader
    _geometries[unit_ii]->setVertexAttribArray(
        2,
        _normalsArrays[unit_ii],
        (_vt_nm != nullptr ? osg::Array::BIND_PER_VERTEX
                           : osg::Array::BIND_PER_PRIMITIVE_SET));

    // if(_colorsArrays[unit_ii].get()->empty() &&
    if((_vt_uv_m != nullptr || _het_uv_m != nullptr) &&
       !_texcoordsArrays[unit_ii].get()->empty())
    {
      program = loadTexturedMesh(_geometries,
                                 _vertexArrays,
                                 _normalsArrays,
                                 _tangentsArrays,
                                 _texcoordsArrays,
                                 m_SmoothFlat_Shading,
                                 _vt_nm,
                                 _vt_tm,
                                 _m_mm,
                                 unit_ii);
    }
    else
    {
      colorsArrays2.push_back(new osg::Vec4Array);

      const std::vector< osg::ref_ptr< osg::Vec4Array > > *colorsArrays_tmp;

      VertexColorMap vt_cm2;
      VertexColorMap *vt_cm_tmp;

      if(_vt_cm != nullptr || _f_cm != nullptr)
      {
        colorsArrays_tmp = &_colorsArrays;
        vt_cm_tmp = _vt_cm;
      }
      else
      {
        colorsArrays2[unit_ii].get()->resize(
            _vertexArrays[unit_ii].get()->size(),
            Helpers::ColorConverter(Color::Wisteria()));

        colorsArrays_tmp = &colorsArrays2;

        vt_cm2 = make_property_map(FEVV::vertex_color, *_g);
        vt_cm_tmp = &vt_cm2;
      }

      program = loadColoredMesh(_geometries,
                                _vertexArrays,
                                _normalsArrays,
                                *colorsArrays_tmp /*_colorsArrays*/,
                                _vt_nm,
                                vt_cm_tmp /*_vt_cm*/,
                                unit_ii);
    }

    _geometries[unit_ii]->getOrCreateStateSet()->setAttributeAndModes(
        program.get(), osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);

    size_t nb_faces = size_of_faces(*_g);

    // attach geometry to geode
    if(nb_faces != 0) // not to be done for 'only_pts' mode
      _geode->addDrawable(_geometries[unit_ii]);

    if(unit_ii == 0) // MUST be done only one time
      if(m_RenderSuperimposedEdges || m_RenderSuperimposedVertices ||
         m_RenderSuperimposedVertices_Big || m_Show_Vertex_Normals ||
         (nb_faces == 0)) // last test for 'only_pts' mode
      {
        // RM: adding predefined basic shaders to display superimposed features
        osg::ref_ptr< osg::Program > superimpProgram = new osg::Program;

        osg::ref_ptr< osg::Shader > superimpVertShader =
            new osg::Shader(osg::Shader::VERTEX);
        if(!superimpVertShader->loadShaderSourceFromFile(
               shadersDirectory + "superimposeVert.glsl"))
          std::cerr << "[MeshLoading] Could not read VERTEX shader from file"
                    << std::endl;
        superimpProgram->addShader(superimpVertShader);

        // uncomment to use GEOMETRY shader
          /*osg::ref_ptr< osg::Shader > superimpGeomShader =
              new osg::Shader(osg::Shader::GEOMETRY);
          if(!superimpGeomShader->loadShaderSourceFromFile(
                 shadersDirectory + "superimposeGeom.glsl"))
            std::cerr << "[MeshLoading] Could not read GEOMETRY shader from file"
                      << std::endl;
          superimpProgram->addShader(superimpGeomShader);*/
        // uncomment to use GEOMETRY shader

        osg::ref_ptr< osg::Shader > superimpFragShader =
            new osg::Shader(osg::Shader::FRAGMENT);
        if(!superimpFragShader->loadShaderSourceFromFile(
               shadersDirectory + "superimposeFrag.glsl"))
          std::cerr << "Could not read FRAGMENT shader from file" << std::endl;
        superimpProgram->addShader(superimpFragShader);

        // RM: loading superimposed features
        // edges - superimpose only
        if(m_RenderSuperimposedEdges &&
           (nb_faces != 0)) // not to be done for 'only_pts' mode
        {
          std::cout << "[MeshLoading] Drawing superimposed edges" << std::endl;

          // RM: disabled these calls, fixing loading time & FPS problems
          //_geometries_edges[unit_ii]->setUseVertexBufferObjects(true);

          _geometries_edges[unit_ii]->setVertexAttribArray(
              0, _vertexArrays_edges[unit_ii], osg::Array::BIND_PER_VERTEX);
          _geometries_edges[unit_ii]->setVertexAttribArray(
              1, _colorsArrays_edges[unit_ii], osg::Array::BIND_PER_VERTEX);
          _geometries_edges[unit_ii]->setVertexAttribArray(
              2, _normalsArrays_edges[unit_ii], osg::Array::BIND_PER_VERTEX);
          _geometries_edges[unit_ii]
              ->getOrCreateStateSet()
              ->setAttributeAndModes(superimpProgram.get(),
                                     osg::StateAttribute::ON |
                                         osg::StateAttribute::OVERRIDE);

          _geode->addDrawable(_geometries_edges[unit_ii]);
        }

        // vertices - superimpose and 'only_pts' mode only
        if(m_RenderSuperimposedVertices || m_RenderSuperimposedVertices_Big ||
           (nb_faces == 0)) // last test for 'only_pts' mode
        {
          std::cout << "[MeshLoading] Drawing superimposed vertices"
                    << std::endl;

          // RM: disabled these calls, fixing loading time & FPS problems
          //_geometries_vertices[unit_ii]->setUseVertexBufferObjects(true);

          _geometries_vertices[unit_ii]->setVertexAttribArray(
              0, _vertexArrays_vertices[unit_ii], osg::Array::BIND_PER_VERTEX);
          _geometries_vertices[unit_ii]->setVertexAttribArray(
              1, _colorsArrays_vertices[unit_ii], osg::Array::BIND_PER_VERTEX);
          _geometries_vertices[unit_ii]->setVertexAttribArray(
              2, _normalsArrays_vertices[unit_ii], osg::Array::BIND_PER_VERTEX);
          _geometries_vertices[unit_ii]
              ->getOrCreateStateSet()
              ->setAttributeAndModes(superimpProgram.get(),
                                     osg::StateAttribute::ON |
                                         osg::StateAttribute::OVERRIDE);

          // RM: if bigger vertices are needed, send a default size of 5 to the
          // shader, which will display them accordingly
          /*if (m_RenderSuperimposedVertices_Big)
          {
            _geometries_vertices[unit_ii]->getOrCreateStateSet()->setMode(GL_PROGRAM_POINT_SIZE,
          osg::StateAttribute::ON);
            _geometries_vertices[unit_ii]->getOrCreateStateSet()->addUniform(new
          osg::Uniform("uniPointSize", 5u));
          }*/

          _geode->addDrawable(_geometries_vertices[unit_ii]);
        }

        // [ normals
        if(m_Show_Vertex_Normals)
        {
          std::cout << "[MeshLoading] Drawing normals"
                    << std::endl;

          // disabled these calls, fixing loading time & FPS problems
          //_geometries_normals[unit_ii]->setUseVertexBufferObjects(true);

          _geometries_normals[unit_ii]->setVertexAttribArray(
              0, _vertexArrays_normals[unit_ii], osg::Array::BIND_PER_VERTEX);
          _geometries_normals[unit_ii]->setVertexAttribArray(
              1, _colorsArrays_normals[unit_ii], osg::Array::BIND_PER_VERTEX);
          _geometries_normals[unit_ii]
              ->getOrCreateStateSet()
              ->setAttributeAndModes(superimpProgram.get(),
                                     osg::StateAttribute::ON |
                                         osg::StateAttribute::OVERRIDE);

          _geode->addDrawable(_geometries_normals[unit_ii]);
        }
        // [ normals
      }

    ++unit_ii;
  } while(unit_ii < _m_mm_size);

  // RM: enable/disable lighting (enabled by default)
  if(m_Lighting)
    _geode->getOrCreateStateSet()->addUniform(
        new osg::Uniform("uniUseLighting", true));
  else
    _geode->getOrCreateStateSet()->addUniform(
        new osg::Uniform("uniUseLighting", false));

  if(!m_SmoothFlat_Shading)
    _geode->getOrCreateStateSet()->addUniform(
        new osg::Uniform("uniUseSmoothShading", false));

  std::cout << "--- Setting uniforms callbacks" << std::endl;

  // osgViewer::View* view = getView(0); // for osgViewer::CompositeViewer
  osgViewer::View *view =
      dynamic_cast< osgViewer::View * >(this); // for osgViewer::Viewer

  // RM: setting callback for camera position
  osg::Uniform *cameraPos =
      new osg::Uniform(osg::Uniform::FLOAT_VEC3, "uniCameraPos");
  cameraPos->setUpdateCallback(new CameraPosCallback(view->getCamera()));
  _geode->getOrCreateStateSet()->addUniform(cameraPos);

  // RM: setting callback for camera target (unneeded for now since view matrix
  // seems to embed it (theorically is a look-at))
  /*osg::Uniform* cameraTarget = new osg::Uniform(osg::Uniform::FLOAT_VEC3,
  "uniCameraTarget"); cameraTarget->setUpdateCallback(new
  CameraTargetCallback(_geode));
  _geode->getOrCreateStateSet()->addUniform(cameraTarget);*/

  // RM: setting callback for model matrix
  osg::Uniform *modelMatrix =
      new osg::Uniform(osg::Uniform::FLOAT_MAT4, "uniModelMatrix");
  modelMatrix->setUpdateCallback(new ModelMatrixCallback());
  _geode->getOrCreateStateSet()->addUniform(modelMatrix);

  // RM: setting callback for view matrix
  osg::Uniform *viewMatrix =
      new osg::Uniform(osg::Uniform::FLOAT_MAT4, "uniViewMatrix");
  viewMatrix->setUpdateCallback(new ViewMatrixCallback(view->getCamera()));
  _geode->getOrCreateStateSet()->addUniform(viewMatrix);

  // RM: setting callback for inverse view matrix
  osg::Uniform *invViewMatrix =
      new osg::Uniform(osg::Uniform::FLOAT_MAT4, "uniInvViewMatrix");
  invViewMatrix->setUpdateCallback(
      new InverseViewMatrixCallback(view->getCamera()));
  _geode->getOrCreateStateSet()->addUniform(invViewMatrix);

  // RM: setting callback for MVP matrix
  osg::Uniform *mvpMatrix =
      new osg::Uniform(osg::Uniform::FLOAT_MAT4, "uniMvpMatrix");
  mvpMatrix->setUpdateCallback(new MVPMatrixCallback(view->getCamera()));
  _geode->getOrCreateStateSet()->addUniform(mvpMatrix);


  // MTO: setting callback for model matrix
  osg::Uniform *g_modelMatrix =
      new osg::Uniform(osg::Uniform::FLOAT_MAT4, "model");
  g_modelMatrix->setUpdateCallback(new ModelMatrixCallback());
  _geode->getOrCreateStateSet()->addUniform(g_modelMatrix);

  // MTO: setting callback for view matrix
  osg::Uniform *g_viewMatrix =
      new osg::Uniform(osg::Uniform::FLOAT_MAT4, "view");
  g_viewMatrix->setUpdateCallback(new ViewMatrixCallback(view->getCamera()));
  _geode->getOrCreateStateSet()->addUniform(g_viewMatrix);

  // MTO: setting callback for projection matrix
  osg::Uniform *g_projectionMatrix =
      new osg::Uniform(osg::Uniform::FLOAT_MAT4, "projection");
  g_projectionMatrix->setUpdateCallback(new ProjectionMatrixCallback(view->getCamera()));
  _geode->getOrCreateStateSet()->addUniform(g_projectionMatrix);

  /*
  RM: creating light(s)
    - Point light: position = x/y/z/1, no need to set direction
    - Directional light: position = x/y/z/0, direction = x/y/z
        (actually we don't really care about directional lights's position
  x/y/z, except if it's to be displayed) (Stressing out: point light's
  direction.w *must be* 1 (actually anything different than 0, but should be 1
  for consistency), directional's *must be* 0)

  Example:
      osg::ref_ptr<osg::Light> new_light = new osg::Light;

      new_light->setPosition(osg::Vec4(0.f, 1.f, 0.f, 1.f)); // The 1.f at the
  end means that this is *not* a directional light
      //new_light->setDirection(osg::Vec3(-1.f, -1.f, 0.f)); // Must be set for
  directional light new_light->setDiffuse(osg::Vec4(1.f, 0.f, 0.f, 1.f)); //
  Light's color (R/G/B), this one will then be red (set to 1.f/1.f/1.f to get a
  white color); the last parameter has no influence
      new_light->setLinearAttenuation(1.f); // Light intensity

      lights.push_back(new_light);
  */

  // RM: sending lights to the shader
  //   This is called for every mesh; ideally, this should be done upstream
  //   *once* with a Uniform Buffer Object
  _geode->getOrCreateStateSet()->addUniform(new osg::Uniform(
      "uniLightCount", static_cast< unsigned int >(lights.size())));

  for(std::size_t lightIndex = 0; lightIndex < lights.size(); ++lightIndex)
  {
    const std::string locationBase =
        "uniLights[" + std::to_string(lightIndex) + "].";

    const std::string posLocation = locationBase + "position";
    const std::string dirLocation = locationBase + "direction";
    const std::string colorLocation = locationBase + "color";
    const std::string energyLocation = locationBase + "energy";

    _geode->getOrCreateStateSet()->addUniform(new osg::Uniform(
        posLocation.c_str(), lights[lightIndex]->getPosition()));
    _geode->getOrCreateStateSet()->addUniform(new osg::Uniform(
        dirLocation.c_str(), lights[lightIndex]->getDirection()));
    _geode->getOrCreateStateSet()->addUniform(new osg::Uniform(
        colorLocation.c_str(), lights[lightIndex]->getDiffuse()));
    _geode->getOrCreateStateSet()->addUniform(new osg::Uniform(
        energyLocation.c_str(), lights[lightIndex]->getLinearAttenuation()));
  }

  osgUtil::Optimizer optimizer; // https://github.com/openscenegraph/OpenSceneGraph/blob/master/include/osgUtil/Optimizer
  /*ALL_OPTIMIZATIONS = FLATTEN_STATIC_TRANSFORMS_DUPLICATING_SHARED_SUBGRAPHS |
                      REMOVE_REDUNDANT_NODES |
                      REMOVE_LOADED_PROXY_NODES |
                      COMBINE_ADJACENT_LODS |
                      SHARE_DUPLICATE_STATE |
                      MERGE_GEODES |
                      MERGE_GEOMETRY |
                      MAKE_FAST_GEOMETRY |
                      CHECK_GEOMETRY |
                      SPATIALIZE_GROUPS |
                      COPY_SHARED_NODES |
                      TRISTRIP_GEOMETRY |
                      OPTIMIZE_TEXTURE_SETTINGS |
                      TEXTURE_ATLAS_BUILDER |
                      STATIC_OBJECT_DETECTION |
                      BUFFER_OBJECT_SETTINGS*/

#if OSG_MIN_VERSION_REQUIRED(3, 6, 0)
  optimizer.optimize(_geode, osgUtil::Optimizer::BUFFER_OBJECT_SETTINGS /*| osgUtil::Optimizer::MERGE_GEODES*/ | osgUtil::Optimizer::MERGE_GEOMETRY /*| osgUtil::Optimizer::MAKE_FAST_GEOMETRY*/);
#else
  optimizer.optimize(_geode, osgUtil::Optimizer::MERGE_GEOMETRY);
#endif
}

template< typename HalfedgeGraph,
          typename VertexNormalMap,
          typename VertexColorMap,
          typename FaceColorMap,
          typename VertexUVMap,
          typename HalfedgeUVMap,
          typename FaceMaterialMap >
void
FEVV::SimpleViewer::internal_loadLegacyMesh(
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
    FaceMaterialMap *_m_mm)
{
  std::cout << "[MeshLoading] Loading using legacy rendering." << std::endl;

  size_t unit_ii = 0;
  do
  {
    _geometries[unit_ii]->setStateSet(NULL);   // NEW
    _geometries[unit_ii]->setColorArray(NULL); // NEW

    _geometries[unit_ii]->setVertexArray(_vertexArrays[unit_ii]);
    if(_vt_nm != nullptr)
    {
      _geometries[unit_ii]->setNormalArray(_normalsArrays[unit_ii],
                                           osg::Array::BIND_PER_VERTEX);
    }
    else
    {
      _geometries[unit_ii]->setNormalArray(_normalsArrays[unit_ii],
                                           osg::Array::BIND_PER_PRIMITIVE_SET);
    }
    if(_vt_cm != nullptr)
    {
      _geometries[unit_ii]->setColorArray(_colorsArrays[unit_ii],
                                          osg::Array::BIND_PER_VERTEX);
    }
    else if(_f_cm != nullptr)
    {
      _geometries[unit_ii]->setColorArray(_colorsArrays[unit_ii],
                                          osg::Array::BIND_PER_PRIMITIVE_SET);
    }
    // TEXTURES
    else if((_het_uv_m != nullptr || _vt_uv_m != nullptr) &&
            (_texture_type == VERTEX_TEXCOORDS2D ||
             _texture_type == HALFEDGE_TEXCOORDS2D))
    {
      auto material = get(*_m_mm, unit_ii);

      osg::ref_ptr< osg::Material > Kdmaterial = new osg::Material;
      Kdmaterial->setDiffuse(osg::Material::FRONT,
                             osg::Vec4(material.diffuse_red_component,
                                       material.diffuse_green_component,
                                       material.diffuse_blue_component,
                                       1.));
      _geometries[unit_ii]->getOrCreateStateSet()->setAttribute(
          Kdmaterial /*, osg::StateAttribute::OVERRIDE*/);

      if(!material.diffuse_texture_filename.empty())
      {
        // load an image by reading a file
        const osg::ref_ptr< osg::Image > image =
            osgDB::readImageFile(material.diffuse_texture_filename);
        if(image)
        {
          // std::cout << "[MeshLoading] Texture file '" << _texture_file << "'
          // loaded." << std::endl;

          osg::ref_ptr< osg::Texture2D > texture = new osg::Texture2D;

          // protect from being optimized away as static state
          texture->setDataVariance(osg::Object::STATIC); // previously DYNAMIC

          // assign the texture to the image we read from file
          texture->setImage(image);

          texture->setWrap(osg::Texture2D::WRAP_S, osg::Texture2D::REPEAT);
          texture->setWrap(osg::Texture2D::WRAP_T, osg::Texture2D::REPEAT);
          texture->setFilter(
              osg::Texture2D::MIN_FILTER,
              osg::Texture2D::LINEAR_MIPMAP_LINEAR); // previously LINEAR but
                                                     // artifacts !
          texture->setFilter(osg::Texture2D::MAG_FILTER,
                             osg::Texture2D::LINEAR);

          // NEW version - without new StateSet - assign texture unit '0' to the
          // texture we just created and enable the texture
          _geometries[unit_ii]->getOrCreateStateSet()->setAttribute(
              new osg::
                  Material /*, osg::StateAttribute::OVERRIDE*/); // IMPORTANT,
                                                                 // pour
                                                                 // suppression
                                                                 // de l'ancien
                                                                 // material
                                                                 // (OSX PB -
                                                                 // SPACE MODE)
          _geometries[unit_ii]
              ->getOrCreateStateSet()
              ->setTextureAttributeAndModes(
                  0 /*unit_ii*/, texture, osg::StateAttribute::ON);

          _geometries[unit_ii]->setTexCoordArray(0 /*unit_ii*/,
                                                 _texcoordsArrays[unit_ii]);
        }
        else
        {
          std::cerr << "-> [SimpleViewer] Couldn't find texture file: "
                    << material.diffuse_texture_filename << std::endl;
        }
      }
    }
    // MATERIAL
    else
    {
#if 1
      // single side painting

      osg::ref_ptr< osg::Material > Kdmaterial = new osg::Material;
      Kdmaterial->setDiffuse(osg::Material::FRONT,
                             Helpers::ColorConverter(Color::Wisteria()));
      _geometries[unit_ii]->getOrCreateStateSet()->setAttribute(
          Kdmaterial /*, osg::StateAttribute::OVERRIDE*/);
#else
      // both side painting

      osg::ref_ptr< osg::LightModel > lightModel = new osg::LightModel;
      // http://www.glprogramming.com/red/chapter05.html
      // "OpenGL reverses the surface normals for back-facing polygons;
      // typically, this means that the surface normals of visible back- and
      // front-facing polygons face the viewer, rather than pointing away. As a
      // result, all polygons are illuminated correctly. However, these
      // additional operations usually make two-sided lighting perform more
      // slowly than the default one-sided lighting."
      lightModel->setTwoSided(true);
      // http://www.opengl.org/sdk/docs/man2/xhtml/glLightModel.xml
      // "GL_SINGLE_COLOR specifies that a single color is generated from the
      // lighting computation for a vertex. GL_SEPARATE_SPECULAR_COLOR
      // specifies that the specular color computation of lighting be stored
      // separately from the remainder of the lighting computation. The specular
      // color is summed into the generated fragment's color after the
      // application of texture mapping (if enabled). The initial value is
      // GL_SINGLE_COLOR."
      lightModel->setColorControl(osg::LightModel::SINGLE_COLOR);
      // as of 2.9.9 this would be the default
      lightModel->setAmbientIntensity(osg::Vec4(0.2f, 0.2f, 0.2f, 1.0f));
      // http://www.glprogramming.com/red/chapter05.html
      // "A local viewpoint tends to yield more realistic results, but since the
      // direction has to be calculated for each vertex, overall performance is
      // decreased with a local viewpoint. By default, an infinite viewpoint is
      // assumed."
      lightModel->setLocalViewer(false);

      geode->getOrCreateStateSet()->setAttributeAndModes(
          lightModel.get()); // now not geode but geometries[unit_ii] !!!
#endif
    }

    size_t nb_faces = size_of_faces(*_g);

    // attach geometry to geode
    if(nb_faces != 0) // not to be done for 'only_pts' mode
      _geode->addDrawable(_geometries[unit_ii]);

    // edges - superimpose only
    if(unit_ii == 0) // MUST be done only one time
      if(m_RenderSuperimposedEdges &&
         (nb_faces != 0)) // not to be done for 'only_pts' mode
      {
        _geometries_edges[unit_ii]->setVertexArray(
            _vertexArrays_edges[unit_ii]);
        _geometries_edges[unit_ii]->setColorArray(_colorsArrays_edges[unit_ii],
                                                  osg::Array::BIND_PER_VERTEX);
        _geode->addDrawable(_geometries_edges[unit_ii]);
      }

    // vertices - superimpose and 'only_pts' mode only
    if(unit_ii == 0) // MUST be done only one time
      if(m_RenderSuperimposedVertices || m_RenderSuperimposedVertices_Big ||
         (nb_faces == 0)) // last test for 'only_pts' mode
      {
        _geometries_vertices[unit_ii]->setVertexArray(
            _vertexArrays_vertices[unit_ii]);
        _geometries_vertices[unit_ii]->setColorArray(
            _colorsArrays_vertices[unit_ii], osg::Array::BIND_PER_VERTEX);
        _geode->addDrawable(_geometries_vertices[unit_ii]);
      }

    // [ normals
    if(unit_ii == 0) // MUST be done only one time
      if(m_Show_Vertex_Normals)
      {
        _geometries_normals[unit_ii]->setVertexArray(
            _vertexArrays_normals[unit_ii]);
        _geometries_normals[unit_ii]->setColorArray(
            _colorsArrays_normals[unit_ii], osg::Array::BIND_PER_VERTEX);
        _geode->addDrawable(_geometries_normals[unit_ii]);
      }
    // ] normals

    unit_ii++;
  } while(unit_ii < _m_mm_size);

  _geode->getOrCreateStateSet()->setMode(
      GL_NORMALIZE, osg::StateAttribute::ON); // IMPORTANT : this ensures
                                              // correct lighting (normals)

  if(m_Lighting)
    _geode->getOrCreateStateSet()->setMode(GL_LIGHTING,
                                           osg::StateAttribute::ON);
  else
    _geode->getOrCreateStateSet()->setMode(GL_LIGHTING,
                                           osg::StateAttribute::OFF);

  if(m_SmoothFlat_Shading)
  {
    osg::ref_ptr< osg::ShadeModel > shadeModel =
        new osg::ShadeModel(osg::ShadeModel::SMOOTH);
    _geode->getOrCreateStateSet()->setAttribute(shadeModel);
  }
  else
  {
    osg::ref_ptr< osg::ShadeModel > shadeModel =
        new osg::ShadeModel(osg::ShadeModel::FLAT);
    _geode->getOrCreateStateSet()->setAttribute(shadeModel);
  }

  osgUtil::Optimizer optimizer; // https://github.com/openscenegraph/OpenSceneGraph/blob/master/include/osgUtil/Optimizer

#if OSG_MIN_VERSION_REQUIRED(3, 6, 0)
  optimizer.optimize(_geode, osgUtil::Optimizer::BUFFER_OBJECT_SETTINGS /*| osgUtil::Optimizer::MERGE_GEODES*/ | osgUtil::Optimizer::MERGE_GEOMETRY /*| osgUtil::Optimizer::MAKE_FAST_GEOMETRY*/);
#else
  optimizer.optimize(_geode, osgUtil::Optimizer::MERGE_GEOMETRY);
#endif
}
