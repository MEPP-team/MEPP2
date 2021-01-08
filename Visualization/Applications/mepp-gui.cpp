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

#include "mepp-gui.h"

#include "FEVV/DataStructures/DataStructures.h"

#include <fstream>
#include <functional>

#include "Visualization/SimpleApplication.h"
#include "Visualization/SimpleWindow.h"
//#include "Visualization/SimpleAdapterVisu.h" // REMOVE ?
#include "Visualization/SimpleViewer.h"

#include "Base/Color.hpp"
#include "Base/Block.h"

// OSG mesh loader
#include "Visualization/Helpers/OSGHelpers.h"

// IO tools
#include "FEVV/Tools/IO/FileUtilities.hpp"

enum OpenWith {
  OPEN_WITH_NONE,
  OPEN_WITH_ALL,
  OPEN_WITH_POLYHEDRON,
  OPEN_WITH_SURFACEMESH,
  OPEN_WITH_LCC,
  OPEN_WITH_OPENMESH,
  OPEN_WITH_AIF,
  OPEN_WITH_CGALPOINTSET,
  OPEN_WITH_PCLPOINTCLOUD
};


int
main(int argc, char **argv)
{
  bool test = false;
  OpenWith mesh_open_with = OPEN_WITH_ALL;
  OpenWith cloud_open_with = OPEN_WITH_ALL;
  std::vector< std::string > mesh_filenames;
  std::vector< std::string > cloud_filenames;
  std::vector< std::string > mesh_extensions = { ".off", ".obj" };
  std::vector< std::string > cloud_extensions = { ".xyz", ".ply", ".pcd" };

  if((argc == 1) || (argv[1] && (strncmp("-psn", argv[1], 4) == 0)))
                    // very specific ".app/GUI" stuff from Apple...
                    // https://discussions.apple.com/thread/1571921?tstart=0
  {
    std::cout << "Usage: " << std::endl
              << "  " << argv[0]
              << " [mesh_filename] [cloud_filename]"
                 " [-test] [-emptytest] [-all|-polyhedron|-surfacemesh|-lcc|-openmesh|-aif|-cgalps|-pclpc]"
              << std::endl
              << "Examples:" << std::endl
              << argv[0] << " casting.off" << std::endl
              << argv[0] << " casting.off casting.xyz" << std::endl
              << argv[0] << " casting.off -test" << std::endl
              << argv[0] << " casting.off -polyhedron casting.xyz" << std::endl
              << argv[0] << " casting.off -polyhedron casting.xyz -pclpc -test" << std::endl;
    // exit(EXIT_FAILURE);
  }
  else if(argc >= 2)
  {
    // parse arguments
    for(int i = 1; i < argc; i++)
    {
      std::string argvi(argv[i]);

      if(argvi == "-test")
      {
        test = true;
      }
      else if(argvi == "-emptytest")
      {
        test = true;
        mesh_open_with = OPEN_WITH_NONE;
        cloud_open_with = OPEN_WITH_NONE;
      }
      else if(argvi == "-polyhedron")
      {
        mesh_open_with = OPEN_WITH_POLYHEDRON;
      }
      else if(argvi == "-surfacemesh")
      {
        mesh_open_with = OPEN_WITH_SURFACEMESH;
      }
      else if(argvi == "-lcc")
      {
        mesh_open_with = OPEN_WITH_LCC;
      }
      else if(argvi == "-openmesh")
      {
        mesh_open_with = OPEN_WITH_OPENMESH;
      }
      else if(argvi == "-aif")
      {
        mesh_open_with = OPEN_WITH_AIF;
      }
      else if(argvi == "-cgalps")
      {
        cloud_open_with = OPEN_WITH_CGALPOINTSET;
      }
      else if(argvi == "-pclpc")
      {
        cloud_open_with = OPEN_WITH_PCLPOINTCLOUD;
      }
      else if(argvi == "-all")
      {
        mesh_open_with = OPEN_WITH_ALL;
        cloud_open_with = OPEN_WITH_ALL;
      }
      else
      {
        // this is a file name
        if(FEVV::FileUtils::has_extension(argvi, mesh_extensions))
          mesh_filenames.push_back(argvi);
        else if(FEVV::FileUtils::has_extension(argvi, cloud_extensions))
          cloud_filenames.push_back(argvi);
      }
    }
  }

// init GUI

#ifdef Q_WS_X11
#if QT_VERSION >= 0x040800
  // Required for multithreaded QGLWidget on Linux/X11,
  // see http://blog.qt.io/blog/2011/06/03/threaded-opengl-in-4-8/
  QApplication::setAttribute(
      Qt::AA_X11InitThreads); // must be called before QApplication app( argc,
                              // argv );
#endif
#endif

  FEVV::Block::begin("init", "Init Visualization.");
  FEVV::SimpleApplication app(argc, argv);
  FEVV::SimpleWindow gui;

#ifndef FEVV_USE_PCL
  #define PCL_VERSION_PRETTY "(no PCL)"
#endif

  gui.setWindowTitle(
      QObject::tr("%1 - %2 - %3 - %4 - Qt (compiled) %5 - Qt (run-time) %6 - "
                  "OSG %7 - CGAL %8 (%9.%10.%11) - PCL %12")
          .arg(MAINWINDOW_TITLE)
          .arg(MEPP_VERSION)
          .arg(ARCHITECTURE)
          .arg(BUILD_TYPE)
          .arg(QT_VERSION_STR)
          .arg(qVersion())
          .arg(osgGetVersion())
          .arg(CGAL_VERSION_STR)
          .arg(CGAL_VERSION_MAJOR)
          .arg(CGAL_VERSION_MINOR)
          .arg(CGAL_VERSION_PATCH)
          .arg(PCL_VERSION_PRETTY));

  gui.init(test);

  gui.loadQtPlugins();

  FEVV::Block::end("init");

#ifdef DEBUG_VISU2
  std::cout << "*** file " << __FILE__ << " line " << __LINE__ << std::endl;
#endif

  // open mesh file(s) provided on command line

#ifdef DEBUG_VISU2
  std::cout << "*** file " << __FILE__ << " line " << __LINE__ << std::endl;
#endif

  if(! mesh_filenames.empty())
  {
#ifdef FEVV_USE_CGAL
    ///// Polyhedron
    if(mesh_open_with == OPEN_WITH_POLYHEDRON || mesh_open_with == OPEN_WITH_ALL)
    {
      FEVV::Block::begin("loading-polyhedron", "Loading Polyhedron mesh.");
      {
        gui.open_SPACE_TIME< FEVV::MeshPolyhedron >(nullptr, mesh_filenames);
      }
      FEVV::Block::end("loading-polyhedron");
    }

    ///// Surface_mesh
    if(mesh_open_with == OPEN_WITH_SURFACEMESH || mesh_open_with == OPEN_WITH_ALL)
    {
      FEVV::Block::begin("loading-surface", "Loading SurfaceMesh mesh.");
      {
        gui.open_SPACE_TIME< FEVV::MeshSurface >(nullptr, mesh_filenames);
      }
      FEVV::Block::end("loading-surface");
      /// [Snippet Displaying Surface SurfaceMesh]
    }

    ///// Linear_cell_complex
    if(mesh_open_with == OPEN_WITH_LCC || mesh_open_with == OPEN_WITH_ALL)
    {
      FEVV::Block::begin("loading-lcc", "Loading LCC mesh.");
      {
        gui.open_SPACE_TIME< FEVV::MeshLCC >(nullptr, mesh_filenames);
      }
      FEVV::Block::end("loading-lcc");
    }
#endif // FEVV_USE_CGAL

#ifdef DEBUG_VISU2
    std::cout << "*** file " << __FILE__ << " line " << __LINE__ << std::endl;
#endif

#ifdef FEVV_USE_OPENMESH
    ///// OpenMesh
    if(mesh_open_with == OPEN_WITH_OPENMESH || mesh_open_with == OPEN_WITH_ALL)
    {
      FEVV::Block::begin("loading-openmesh", "Loading OpenMesh mesh.");
      {
        gui.open_SPACE_TIME< FEVV::MeshOpenMesh >(nullptr, mesh_filenames);
      }
      FEVV::Block::end("loading-openmesh");
    }
#endif // FEVV_USE_OPENMESH

#ifdef DEBUG_VISU2
    std::cout << "*** file " << __FILE__ << " line " << __LINE__ << std::endl;
#endif

#ifdef FEVV_USE_AIF
    ///// AIF
    if(mesh_open_with == OPEN_WITH_AIF || mesh_open_with == OPEN_WITH_ALL)
    {
      FEVV::Block::begin("loading-aif", "Loading AIF mesh.");
      {
        gui.open_SPACE_TIME< FEVV::MeshAIF >(nullptr, mesh_filenames);
      }
      FEVV::Block::end("loading-aif");
    }
#endif // FEVV_USE_AIF
  }

  // open point cloud file(s) provided on command line

  if(! cloud_filenames.empty())
  {
#ifdef FEVV_USE_CGAL
    // CGALPointSet
    if(cloud_open_with == OPEN_WITH_CGALPOINTSET || cloud_open_with == OPEN_WITH_ALL)
    {
      FEVV::Block::begin("loading-CGALPointSet", "Loading CGALPointSet point cloud.");
      {
        gui.open_SPACE_TIME< FEVV::CGALPointSet >(nullptr, cloud_filenames);
      }
      FEVV::Block::end("loading-CGALPointSet");
    }
#endif // FEVV_USE_CGAL

#ifdef FEVV_USE_PCL
    // PCLPointCloud
    if(cloud_open_with == OPEN_WITH_PCLPOINTCLOUD || cloud_open_with == OPEN_WITH_ALL)
    {
      FEVV::Block::begin("loading-PCLPointCloud", "Loading PCLPointCloud point cloud.");
      {
        gui.open_SPACE_TIME< FEVV::PCLPointCloud >(nullptr, cloud_filenames);
      }
      FEVV::Block::end("loading-PCLPointCloud");
    }
#endif // FEVV_USE_PCL
  }

#ifdef DEBUG_VISU2
  std::cout << "*** file " << __FILE__ << " line " << __LINE__ << std::endl;
#endif

  // run GUI

  gui.show();

#ifdef DEBUG_VISU2
  std::cout << "*** file " << __FILE__ << " line " << __LINE__ << std::endl;
#endif

  int ret = app.exec();

#ifdef DEBUG_VISU2
  std::cout << "*** file " << __FILE__ << " line " << __LINE__ << std::endl;
#endif

  return ret;
}
