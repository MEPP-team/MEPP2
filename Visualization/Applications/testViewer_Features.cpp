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
#if(_MSC_VER >= 1400)
#ifndef _SCL_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif
#endif

#include "FEVV/DataStructures/DataStructures.h"

#include "Visualization/SimpleApplication.h"
#include "Visualization/SimpleWindow.h"
#include "Visualization/SimpleAdapterVisu.h"
#include "Visualization/SimpleViewer.h"

#include "Base/Color.hpp"
#include "Base/Block.h"
#include "Base/GradientColorMap.h"

#include "Visualization/Helpers/OSGHelpers.h"

#include "FEVV/Filters/Generic/calculate_face_normals.hpp"
#include "FEVV/Filters/Generic/translation.hpp"
#include "FEVV/Filters/Generic/scaling.hpp"
#include "FEVV/Filters/uvs.h"

#include "Visualization/PluginFilters/TranslationPlugin.h"

#include <functional>
#include <fstream>

using namespace FEVV;
using namespace FEVV::Filters;
using Mesh = MeshOpenMesh;

// ELO: why do not use 'num_vertices(g)'
//     from
//     http://www.boost.org/doc/libs/1_60_0/libs/graph/doc/graph_concepts.html ?
template< typename HalfedgeGraph >
unsigned int
calculate_vertex_number(const HalfedgeGraph &_g)
{
  using GraphTraits = boost::graph_traits< HalfedgeGraph >;
  using vertex_iterator = typename GraphTraits::vertex_iterator;

  unsigned int size = 0;

  for(vertex_iterator v_it = vertices(_g).first; v_it != vertices(_g).second;
      ++v_it)
  {
    ++size;
  }

  return size;
}

template< typename HalfedgeGraph >
unsigned int
calculate_face_number(const HalfedgeGraph &_g)
{
  using GraphTraits = boost::graph_traits< HalfedgeGraph >;
  using face_iterator = typename GraphTraits::face_iterator;

  unsigned int size = 0;

  for(face_iterator f_it = faces(_g).first; f_it != faces(_g).second; ++f_it)
  {
    ++size;
  }

  return size;
}

template< typename HalfedgeGraph,
          typename ColorMap,
          typename CustomColorMap = FEVV::GradientColorMap< unsigned int > >
void
calculate_gradient_colormap(const HalfedgeGraph &_g,
                            ColorMap _cm,
                            const CustomColorMap &_ccm)
{
  using GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph >;
  using GraphTraits = boost::graph_traits< HalfedgeGraph >;
  using vertex_iterator = typename GraphTraits::vertex_iterator;
  using Vector = typename GeometryTraits::Vector;

  unsigned int size = 0;
  Color c;

  for(vertex_iterator v_it = vertices(_g).first; v_it != vertices(_g).second;
      ++v_it)
  {
    c = _ccm(size);
    _cm[*v_it] = Vector(c.red(), c.green(), c.blue());
    ++size;
  }
}

template< typename HalfedgeGraph,
          typename ColorMap,
          typename CustomColorMap = FEVV::GradientColorMap< unsigned int > >
void
calculate_gradient_facecolormap(const HalfedgeGraph &_g,
                                ColorMap _cm,
                                const CustomColorMap &_ccm)
{
  using GraphTraits = boost::graph_traits< HalfedgeGraph >;
  using face_iterator = typename GraphTraits::face_iterator;
  using GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph >;
  using Vector = typename GeometryTraits::Vector;

  unsigned int size = 0;
  Color c;

  for(face_iterator f_it = faces(_g).first; f_it != faces(_g).second; ++f_it)
  {
    c = _ccm(size);
    _cm[*f_it] = Vector(c.red(), c.green(), c.blue());
    ++size;
  }
}

Mesh *
loadMesh(std::string &_in)
{
  Mesh *mesh = new Mesh;
  try
  {
    OpenMesh::IO::Options options = OpenMesh::IO::Options::Default;
    bool bClear = false;
    OpenMesh::IO::read_mesh(*mesh, _in, options, bClear);
  }
  catch(const std::length_error &le)
  {
    std::cerr << "[OpenMesh] Exception caught while reading input file: "
              << le.what() << std::endl;
    BOOST_ASSERT_MSG(false,
                     "[OpenMesh] Exception caught while reading input file.");
  }

  return mesh;
}

Mesh *
loadMesh_Tex(std::string &_in, int type = VERTEX_TEXCOORDS2D)
{
  Mesh *mesh = new Mesh;
  OpenMesh::IO::Options options;

  if(type == VERTEX_TEXCOORDS2D) // vertex_texcoords2D
  {
    mesh->request_vertex_texcoords2D();
    options += OpenMesh::IO::Options::VertexTexCoord;
  }
  else if(type == HALFEDGE_TEXCOORDS2D) // halfedge_texcoords2D
  {
    mesh->request_halfedge_texcoords2D();
    options += OpenMesh::IO::Options::FaceTexCoord;
  }

  try
  {
    bool bClear = false;
    OpenMesh::IO::read_mesh(*mesh, _in, options, bClear);
  }
  catch(const std::length_error &le)
  {
    std::cerr << "[OpenMesh] Exception caught while reading input file: "
              << le.what() << std::endl;
    BOOST_ASSERT_MSG(false,
                     "[OpenMesh] Exception caught while reading input file.");
  }

  if(mesh->has_vertex_texcoords2D())
  {
    std::cout << "Mesh provides vertex texture coordinates.\n";
  }
  if(mesh->has_halfedge_texcoords2D())
  {
    std::cout << "Mesh provides halfedge texture coordinates.\n";
  }

  // mesh.release_vertex_texcoords2D();
  // mesh.release_halfedge_texcoords2D();

  return mesh;
}

int
main(int argc, char **argv)
{
  if((argc < 2) ||
     (argv[1] &&
      (strncmp("-psn", argv[1], 4) ==
       0))) // very specific ".app/GUI" stuff from Apple... No comment... -
            // https://discussions.apple.com/thread/1571921?tstart=0
  {
    std::cout << "Usage: ./testViewer_Features filename.off [test]" << std::endl
              << "- filename being an off file." << std::endl
              << "- test." << std::endl;
    exit(EXIT_FAILURE);
  }

  bool test = false;
  bool load = true;

  if(argc == 3 && strcmp(argv[2], "test") == 0)
  {
    test = true;
  }
  else if(argc == 3 && strcmp(argv[2], "emptytest") == 0)
  {
    test = true;
    load = false;
  }

#ifdef Q_WS_X11
#if QT_VERSION >= 0x040800
  // Required for multithreaded QGLWidget on Linux/X11,
  // see http://blog.qt.io/blog/2011/06/03/threaded-opengl-in-4-8/
  QApplication::setAttribute(
      Qt::AA_X11InitThreads); // must be called before QApplication app( argc,
                              // argv );
#endif
#endif

  Block::begin("loading", "Loading Visualization.");
  SimpleApplication app(argc, argv);

  SimpleWindow *gui = new SimpleWindow();

  gui->init(test);
  gui->show();

  SimpleAdapterVisu< Mesh > *adapter = nullptr;
  SimpleViewer< Mesh > *viewer;
  if(load)
  {
    adapter = new SimpleAdapterVisu< Mesh >();
    adapter->setWindowTitle(QObject::tr("<OpenMesh>"));
    adapter->setWindowIcon(QIcon(":/logo/resources/MEPP.png"));
    viewer = new SimpleViewer< Mesh >();
    adapter->attach(viewer);
    adapter->init();
    viewer->init();
    gui->attach(adapter, USE_MDI);

    TranslationPlugin *translation = new TranslationPlugin;
    translation->init(0.0, 0.0, 0.0);
    gui->attachPlugin(translation);
  }

  if(USE_MDI)
  {
    if(gui->getMdiArea())
    {
      gui->getMdiArea()->setActivationOrder(QMdiArea::CreationOrder);
      gui->getMdiArea()->tileSubWindows();
    }
  }

  gui->loadQtPlugins();

  // QCoreApplication::processEvents(); ///<! draw window before huge
  // computation

  Block::end("loading");

  if(load)
  {
    Block::begin("loading-openmesh", "Loading OpenMesh mesh.");
    {
      using VectorOpenMesh = Mesh::Normal;

      using VertexMapOpenMesh =
          boost::property_map< Mesh, boost::vertex_index_t >::const_type;
      using FaceMapOpenMesh =
          boost::property_map< Mesh, boost::face_index_t >::const_type;
      using HalfedgeMapOpenMesh =
          boost::property_map< Mesh, boost::halfedge_index_t >::const_type;

      using FaceNormalMapOpenMesh =
          boost::vector_property_map< VectorOpenMesh, FaceMapOpenMesh >;
      using VertexColorMapOpenMesh =
          boost::vector_property_map< VectorOpenMesh, VertexMapOpenMesh >;
      using VertexUVMapOpenMesh =
          boost::vector_property_map< VectorOpenMesh, VertexMapOpenMesh >;
      using HalfedgeUVMapOpenMesh =
          boost::vector_property_map< VectorOpenMesh, HalfedgeMapOpenMesh >;

      // normal test
      std::string in(argv[1]);

      FEVV::PMapsContainer *p_openmesh_pmaps_bag =
          new FEVV::PMapsContainer; // already destroy by the viewer destructor
                                    // // NOT USED, just here for
                                    // compatibility...

      QApplication::setOverrideCursor(Qt::WaitCursor);
      Mesh *mesh = loadMesh(in);
      QApplication::restoreOverrideCursor();

      int tex_type = NO_TEXCOORDS;
      std::string in_Tex("");

      // manual tests - please don't delete these lines
      // std::string in("C:\\data\\textures\\usb\\usb_vt.ply"); int tex_type =
      // VERTEX_TEXCOORDS2D; Mesh* mesh = loadMesh_Tex(in, tex_type); std::string
      // in_Tex("C:\\data\\textures\\usb\\usb.jpg"); std::string
      // in("C:\\data\\textures\\usb\\usb_vt.obj"); int tex_type =
      // VERTEX_TEXCOORDS2D; Mesh* mesh = loadMesh_Tex(in, tex_type); std::string
      // in_Tex("C:\\data\\textures\\usb\\usb.jpg");

      // std::string in("C:\\data\\textures\\kip\\kip.obj"); int tex_type =
      // HALFEDGE_TEXCOORDS2D; Mesh* mesh = loadMesh_Tex(in, tex_type);
      // std::string in_Tex("C:\\data\\textures\\kip\\Kip.png"); std::string
      // in("C:\\data\\textures\\RensoKuster\\RensoKuster.obj"); int tex_type =
      // HALFEDGE_TEXCOORDS2D; Mesh* mesh = loadMesh_Tex(in, tex_type);
      // std::string
      // in_Tex("C:\\data\\textures\\RensoKuster\\RensoKuster_0.jpg");
      // manual tests - please don't delete these lines

      // point map
      auto pm = get(boost::vertex_point, *mesh);
      VertexMapOpenMesh vm = get(boost::vertex_index, *mesh);
      FaceMapOpenMesh fm = get(boost::face_index, *mesh);
      HalfedgeMapOpenMesh hem = get(boost::halfedge_index, *mesh);

      FaceNormalMapOpenMesh f_nm(fm);
      calculate_face_normals(*mesh, pm, f_nm);

      VertexUVMapOpenMesh vt_uv(vm);
      if(tex_type == VERTEX_TEXCOORDS2D)
      {
        set_vertex_uv(*mesh, vt_uv);
        // old specific OpenMesh version (to delete later...)
        /*{
            MeshOpenMesh::FaceIter f_it, f_end(mOpenMesh.faces_end());
            for (f_it = mOpenMesh.faces_begin(); f_it != f_end; ++f_it)
            {
                for (MeshOpenMesh::CFVIter fv_it = mOpenMesh.cfv_iter(*f_it);
        fv_it.is_valid(); ++fv_it)
                {
                    Vector uv(mOpenMesh.texcoord2D(*fv_it)[0],
        mOpenMesh.texcoord2D(*fv_it)[1], 0); vt_uv[*fv_it] = uv;
                }
            }
        }*/
      }

      HalfedgeUVMapOpenMesh he_uv(hem);
      if(tex_type == HALFEDGE_TEXCOORDS2D)
      {
        set_halfedge_uv(*mesh, he_uv);
        // old specific OpenMesh version (to delete later...)
        /*{
            MeshOpenMesh::FaceIter f_it, f_end(mOpenMesh.faces_end());
            for (f_it = mOpenMesh.faces_begin(); f_it != f_end; ++f_it)
            {
                for (MeshOpenMesh::CFHIter fh_it = mOpenMesh.cfh_iter(*f_it);
        fh_it.is_valid(); ++fh_it)
                {
                    Vector uv(mOpenMesh.texcoord2D(*fh_it)[0],
        mOpenMesh.texcoord2D(*fh_it)[1], 0); het_uv[*fh_it] = uv;
                }
            }
        }*/
      }

      VertexColorMapOpenMesh cm(vm);
      unsigned int nbVertex = calculate_vertex_number(*mesh);
      GradientColorMap< unsigned int > gradient(
          0, nbVertex, Color::White(), Color::Black());
      calculate_gradient_colormap(*mesh, cm, gradient);

      // adapter->attachMesh( mesh );

      if(tex_type == VERTEX_TEXCOORDS2D)
      {
        viewer->drawMesh< decltype(pm),
                          SimpleViewer< Mesh >::DefaultVertexNormalMap,
                          FaceNormalMapOpenMesh,
                          VertexColorMapOpenMesh,
                          SimpleViewer< Mesh >::DefaultFaceColorMap,
                          VertexUVMapOpenMesh,
                          SimpleViewer< Mesh >::DefaultHalfedgeUVMap >(
            mesh,
            p_openmesh_pmaps_bag,
            &pm,
            nullptr,
            &f_nm,
            &cm,
            nullptr,
            &vt_uv,
            nullptr,
            tex_type,
            in_Tex);
      }
      else if(tex_type == HALFEDGE_TEXCOORDS2D)
      {
        viewer->drawMesh< decltype(pm),
                          SimpleViewer< Mesh >::DefaultVertexNormalMap,
                          FaceNormalMapOpenMesh,
                          VertexColorMapOpenMesh,
                          SimpleViewer< Mesh >::DefaultFaceColorMap,
                          SimpleViewer< Mesh >::DefaultVertexUVMap,
                          HalfedgeUVMapOpenMesh >(mesh,
                                                  p_openmesh_pmaps_bag,
                                                  &pm,
                                                  nullptr,
                                                  &f_nm,
                                                  &cm,
                                                  nullptr,
                                                  nullptr,
                                                  &he_uv,
                                                  tex_type,
                                                  in_Tex);
      }
      else
      {
        viewer->drawMesh< decltype(pm),
                          SimpleViewer< Mesh >::DefaultVertexNormalMap,
                          FaceNormalMapOpenMesh,
                          VertexColorMapOpenMesh >(
            mesh, p_openmesh_pmaps_bag, &pm, nullptr, &f_nm, &cm);
      }

      viewer->centerMesh(mesh);
    }
    Block::end("loading-openmesh");

    Block::begin("loading-openmesh-2", "Loading second mesh");
    {
      std::string in(argv[1]);

      FEVV::PMapsContainer *p_openmesh_pmaps_bag =
          new FEVV::PMapsContainer; // already destroy by the viewer destructor
                                    // // NOT USED, just here for
                                    // compatibility...

      QApplication::setOverrideCursor(Qt::WaitCursor);
      Mesh *mesh2 = loadMesh(in);
      QApplication::restoreOverrideCursor();

      auto pm2 = get(boost::vertex_point, *mesh2);

      translate(*mesh2, pm2, -5, 0, 0);

      viewer->drawMesh(mesh2, p_openmesh_pmaps_bag, &pm2);

      viewer->centerMesh(mesh2);
    }
    Block::end("loading-openmesh-2");
  }

  int ret = app.exec();

  // delete viewer; // already destroy by the adapter destructor
  if(!USE_MDI)
    if(adapter != nullptr)
      delete adapter;

  delete gui;

  return ret;
}
