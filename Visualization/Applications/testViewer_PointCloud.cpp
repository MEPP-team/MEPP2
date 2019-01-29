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

#include <fstream>
#include <functional>

#include "Visualization/SimpleApplication.h"
#include "Visualization/SimpleWindow.h"
#include "Visualization/SimpleAdapterVisu.h"
#include "Visualization/SimpleViewer.h"

#include "Base/Color.hpp"
#include "Base/Block.h"

// OSG mesh loader
#include "Visualization/Helpers/OSGHelpers.h"

#if defined _MSC_VER
#pragma warning(disable : 4996) // MT
#endif

#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>

using namespace FEVV;

/// @todo For the moment, we need a model of HalfedgeGraph concept to build GUI.
/// This should be removed soon.
using UnusedMesh = MeshOpenMesh;


int
main(int argc, char **argv)
{
  if((argc < 2) ||
     (argv[1] &&
      (strncmp("-psn", argv[1], 4) ==
       0))) // very specific ".app/GUI" stuff from Apple... No comment... -
            // https://discussions.apple.com/thread/1571921?tstart=0
  {
    std::cout << "Usage: ./testViewer_PointCloud filename.pcd [withColor = "
                 "false] [test]"
              << std::endl
              << "- filename being an pcd file." << std::endl
              << "- withColor being a boolean (0 or 1)." << std::endl
              << "- test." << std::endl;
    exit(EXIT_FAILURE);
  }

  bool test = false;
  bool load = true;

  if(argc == 4 && strcmp(argv[3], "test") == 0)
  {
    test = true;
  }
  else if(argc == 4 && strcmp(argv[3], "emptytest") == 0)
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
  SimpleWindow *gui;

  gui = new SimpleWindow();

  gui->init(test);
  gui->show();

  SimpleAdapterVisu< UnusedMesh > *adapter = nullptr;
  SimpleViewer< UnusedMesh > *viewer;
  if(load)
  {
    adapter = new SimpleAdapterVisu< UnusedMesh >();
    adapter->setWindowTitle(QObject::tr("<PointCloud>"));
    adapter->setWindowIcon(QIcon(":/logo/resources/MEPP.png"));
    viewer = new SimpleViewer< UnusedMesh >();
    adapter->attach(viewer);
    adapter->init();
    viewer->init();
    gui->attach(adapter, USE_MDI);
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
    Block::begin("loading-pcl", "Loading PCL file.");
    {
      std::string in(argv[1]);

      bool with_color = false;
      if(argc >= 3)
      {
        if((strcmp(argv[2], "true") == 0) || (strcmp(argv[2], "True") == 0) ||
           (strcmp(argv[2], "1") == 0))
        {
          with_color = true;
        }
      }

      pcl::PointCloud< pcl::PointXYZRGB >::Ptr cloud(
          new pcl::PointCloud< pcl::PointXYZRGB >);

      if(pcl::io::loadPCDFile< pcl::PointXYZRGB >(in, *cloud) == -1)
      {
        PCL_ERROR("Couldn't read file test_pcd.pcd \n");
        // delete viewer; // already destroy by the adapter destructor
        if(!USE_MDI)
          if(adapter != nullptr)
            delete adapter;
        delete gui;
        return (-1);
      }

      viewer->addGroup(viewer->createMesh(cloud->points, with_color));
    }
    Block::end("loading-pcl");
  }

  int ret = app.exec();

  /////// Cleaning stuff
  // delete viewer; // already destroy by the adapter destructor
  if(!USE_MDI)
    if(adapter != nullptr)
      delete adapter;

  delete gui;
  /////// Cleaning stuff

  return ret;
}
