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

// Qt
#include <QDesktopWidget>

enum Open_with {
  OPEN_WITH_NONE,
  OPEN_WITH_ALL,
  OPEN_WITH_POLYHEDRON,
  OPEN_WITH_SURFACEMESH,
  OPEN_WITH_LCC,
  OPEN_WITH_OPENMESH,
  OPEN_WITH_AIF
};
Open_with open_with = OPEN_WITH_ALL;


int
main(int argc, char **argv)
{
  bool test = false;
  open_with = OPEN_WITH_NONE;

  if((argc == 1) ||
     (argv[1] &&
      (strncmp("-psn", argv[1], 4) ==
       0))) // very specific ".app/GUI" stuff from Apple... No comment... -
            // https://discussions.apple.com/thread/1571921?tstart=0
  {
    std::cout << "Usage: " << argv[0]
              << " filename.off [test] "
                 "[-all][-polyhedron][-surfacemesh][-lcc][-openmesh][-aif]"
              << std::endl
              << "- filename being an off file." << std::endl;
    // exit(EXIT_FAILURE);
  }
  else if(argc == 2)
  {
    open_with = OPEN_WITH_ALL;
  }
  else if(argc >= 3)
  {
    // parse arguments
    for(int i = 2; i < argc; i++)
    {
      std::string argvi(argv[i]);

      if(argvi == "test")
      {
        test = true;
        open_with = OPEN_WITH_ALL;
      }
      else if(argvi == "emptytest")
      {
        test = true;
        open_with = OPEN_WITH_NONE;
      }
      else if(argvi == "-polyhedron")
      {
        open_with = OPEN_WITH_POLYHEDRON;
      }
      else if(argvi == "-surfacemesh")
      {
        open_with = OPEN_WITH_SURFACEMESH;
      }
      else if(argvi == "-lcc")
      {
        open_with = OPEN_WITH_LCC;
      }
      else if(argvi == "-openmesh")
      {
        open_with = OPEN_WITH_OPENMESH;
      }
      else if(argvi == "-aif")
      {
        open_with = OPEN_WITH_AIF;
      }
      else if(argvi == "-all")
      {
        open_with = OPEN_WITH_ALL;
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

  gui.setWindowTitle(
      QObject::tr("%1 - %2 - %3 - %4 - Qt (compiled) %5 - Qt (run-time) %6 - "
                  "OSG %7 - CGAL %8 (%9.%10.%11)")
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
          .arg(CGAL_VERSION_PATCH));

  gui.init(test);

  // gui->resize(QDesktopWidget().availableGeometry().size() * 0.7);
  QRect screen_size = QDesktopWidget().availableGeometry();
  int win_w = screen_size.width() * 0.9;
  int win_h = screen_size.height() * 0.8;
  int pos_x = (screen_size.width() - win_w) / 2;
  int pos_y = (screen_size.height() - win_h) / 2;
  gui.move(pos_x, pos_y);
  gui.resize(win_w, win_h);

  gui.loadQtPlugins();

  FEVV::Block::end("init");

#ifdef DEBUG_VISU2
  std::cout << "*** file " << __FILE__ << " line " << __LINE__ << std::endl;
#endif

  // open mesh file provided on command line

  std::vector< std::string > mesh_filenames;
  if(open_with != OPEN_WITH_NONE)
    mesh_filenames.push_back(std::string(argv[1]));

#ifdef DEBUG_VISU2
  std::cout << "*** file " << __FILE__ << " line " << __LINE__ << std::endl;
#endif

#ifdef FEVV_USE_CGAL
  ///// Polyhedron
  if(open_with == OPEN_WITH_POLYHEDRON || open_with == OPEN_WITH_ALL)
  {
    FEVV::Block::begin("loading-polyhedron", "Loading Polyhedron mesh.");
    {
      gui.open_SPACE_TIME< FEVV::MeshPolyhedron >(nullptr, mesh_filenames);
    }
    FEVV::Block::end("loading-polyhedron");
  }

  ///// Surface_mesh
  if(open_with == OPEN_WITH_SURFACEMESH || open_with == OPEN_WITH_ALL)
  {
    FEVV::Block::begin("loading-surface", "Loading SurfaceMesh mesh.");
    {
      gui.open_SPACE_TIME< FEVV::MeshSurface >(nullptr, mesh_filenames);
    }
    FEVV::Block::end("loading-surface");
    /// [Snippet Displaying Surface SurfaceMesh]
  }

  ///// Linear_cell_complex
  if(open_with == OPEN_WITH_LCC || open_with == OPEN_WITH_ALL)
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
  if(open_with == OPEN_WITH_OPENMESH || open_with == OPEN_WITH_ALL)
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
  if(open_with == OPEN_WITH_AIF || open_with == OPEN_WITH_ALL)
  {
    FEVV::Block::begin("loading-aif", "Loading AIF mesh.");
    {
      gui.open_SPACE_TIME< FEVV::MeshAIF >(nullptr, mesh_filenames);
    }
    FEVV::Block::end("loading-aif");
  }
#endif // FEVV_USE_AIF

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
