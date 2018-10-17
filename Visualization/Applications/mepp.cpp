#if(_MSC_VER >= 1400)
#ifndef _SCL_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif
#endif

#include "mepp.h"

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

#ifdef Q_WS_X11
#if QT_VERSION >= 0x040800
  // Required for multithreaded QGLWidget on Linux/X11,
  // see http://blog.qt.io/blog/2011/06/03/threaded-opengl-in-4-8/
  QApplication::setAttribute(
      Qt::AA_X11InitThreads); // must be called before QApplication app( argc,
                              // argv );
#endif
#endif

  FEVV::Block::begin("loading", "Loading Visualization.");
  FEVV::SimpleApplication app(argc, argv);
  FEVV::SimpleWindow *gui;

  gui = new FEVV::SimpleWindow();
  gui->setWindowTitle(
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

  gui->init(test);

  // gui->resize(QDesktopWidget().availableGeometry().size() * 0.7);
  QRect screen_size = QDesktopWidget().availableGeometry();
  gui->resize(screen_size.width() * 0.9, screen_size.height() * 0.8);

  gui->show();

///////////////////////////////////////////////////////////////////////////////
/// Creation of adapters (reverse order because of MDI view)
///////////////////////////////////////////////////////////////////////////////
#ifdef FEVV_USE_AIF
  ///// AIF
  FEVV::SimpleAdapterVisu< FEVV::MeshAIF > *adapter_aif = nullptr;
  FEVV::SimpleViewer< FEVV::MeshAIF > *viewer_aif;
  if(open_with == OPEN_WITH_AIF || open_with == OPEN_WITH_ALL)
  {
    adapter_aif = new FEVV::SimpleAdapterVisu< FEVV::MeshAIF >();
    // adapter_aif->setWindowTitle(QObject::tr("AIFMesh"));
    adapter_aif->setWindowTitle(
        QObject::tr("<AIFMesh - aid: %1>").arg((qlonglong)adapter_aif, 0, 16));
    adapter_aif->setWindowIcon(QIcon(":/logo/resources/MEPP.png"));
    viewer_aif = new FEVV::SimpleViewer< FEVV::MeshAIF >();
    adapter_aif->attach(viewer_aif);
    adapter_aif->init();
    viewer_aif->init();
    gui->attach(adapter_aif, USE_MDI);
  }
#endif

#ifdef FEVV_USE_OPENMESH
  ///// OpenMesh
  FEVV::SimpleAdapterVisu< FEVV::MeshOpenMesh > *adapter_openmesh = nullptr;
  FEVV::SimpleViewer< FEVV::MeshOpenMesh > *viewer_openmesh;
  if(open_with == OPEN_WITH_OPENMESH || open_with == OPEN_WITH_ALL)
  {
    adapter_openmesh = new FEVV::SimpleAdapterVisu< FEVV::MeshOpenMesh >();
    // adapter_openmesh->setWindowTitle(QObject::tr("OpenMesh"));
    adapter_openmesh->setWindowTitle(
        QObject::tr("<OpenMesh - aid: %1>")
            .arg((qlonglong)adapter_openmesh, 0, 16));
    adapter_openmesh->setWindowIcon(QIcon(":/logo/resources/MEPP.png"));
    viewer_openmesh = new FEVV::SimpleViewer< FEVV::MeshOpenMesh >();
    adapter_openmesh->attach(viewer_openmesh);
    adapter_openmesh->init();
    viewer_openmesh->init();
    gui->attach(adapter_openmesh, USE_MDI);
  }
#endif

#ifdef FEVV_USE_CGAL
  ///// Linear_cell_complex
  FEVV::SimpleAdapterVisu< FEVV::MeshLCC > *adapter_lcc = nullptr;
  FEVV::SimpleViewer< FEVV::MeshLCC > *viewer_lcc;
  if(open_with == OPEN_WITH_LCC || open_with == OPEN_WITH_ALL)
  {
    adapter_lcc = new FEVV::SimpleAdapterVisu< FEVV::MeshLCC >();
    // adapter_lcc->setWindowTitle(QObject::tr("Linear_cell_complex"));
    adapter_lcc->setWindowTitle(QObject::tr("<Linear_cell_complex - aid: %1>")
                                    .arg((qlonglong)adapter_lcc, 0, 16));
    adapter_lcc->setWindowIcon(QIcon(":/logo/resources/MEPP.png"));
    viewer_lcc = new FEVV::SimpleViewer< FEVV::MeshLCC >();
    adapter_lcc->attach(viewer_lcc);
    adapter_lcc->init();
    viewer_lcc->init();
    gui->attach(adapter_lcc, USE_MDI);
  }

  ///// Surface_mesh
  FEVV::SimpleAdapterVisu< FEVV::MeshSurface > *adapter_surface = nullptr;
  FEVV::SimpleViewer< FEVV::MeshSurface > *viewer_surface;
  if(open_with == OPEN_WITH_SURFACEMESH || open_with == OPEN_WITH_ALL)
  {
    adapter_surface = new FEVV::SimpleAdapterVisu< FEVV::MeshSurface >();
    // adapter_surface->setWindowTitle(QObject::tr("Surface_mesh"));
    adapter_surface->setWindowTitle(
        QObject::tr("<Surface_mesh - aid: %1>")
            .arg((qlonglong)adapter_surface, 0, 16));
    adapter_surface->setWindowIcon(QIcon(":/logo/resources/MEPP.png"));
    viewer_surface = new FEVV::SimpleViewer< FEVV::MeshSurface >();
    adapter_surface->attach(viewer_surface);
    adapter_surface->init();
    viewer_surface->init();
    gui->attach(adapter_surface, USE_MDI);
  }

  ///// Polyhedron
  FEVV::SimpleAdapterVisu< FEVV::MeshPolyhedron > *adapter_polyhedron = nullptr;
  FEVV::SimpleViewer< FEVV::MeshPolyhedron > *viewer_polyhedron;
  if(open_with == OPEN_WITH_POLYHEDRON || open_with == OPEN_WITH_ALL)
  {
    adapter_polyhedron = new FEVV::SimpleAdapterVisu< FEVV::MeshPolyhedron >();
    // adapter_polyhedron->setWindowTitle(QObject::tr("Polyhedron_3"));
    adapter_polyhedron->setWindowTitle(
        QObject::tr("<Polyhedron_3 - aid: %1>")
            .arg((qlonglong)adapter_polyhedron, 0, 16));
    adapter_polyhedron->setWindowIcon(QIcon(":/logo/resources/MEPP.png"));
    viewer_polyhedron = new FEVV::SimpleViewer< FEVV::MeshPolyhedron >();
    adapter_polyhedron->attach(viewer_polyhedron);
    adapter_polyhedron->init();
    viewer_polyhedron->init();
    gui->attach(adapter_polyhedron, USE_MDI);
  }
#endif
  ///////////////////////////////////////////////////////////////////////////////

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

  FEVV::Block::end("loading");

  std::string mesh_filename;
  if(open_with != OPEN_WITH_NONE)
    mesh_filename = std::string(argv[1]);

#ifdef FEVV_USE_CGAL
  ///// Polyhedron
  FEVV::MeshPolyhedron *m_polyhedron = nullptr;
  if(open_with == OPEN_WITH_POLYHEDRON || open_with == OPEN_WITH_ALL)
  {
    FEVV::PMapsContainer *p_polyhedron_pmaps_bag =
        new FEVV::PMapsContainer; // already destroy by the viewer destructor
    FEVV::Block::begin("loading-polyhedron", "Loading Polyhedron mesh.");
    {
      QApplication::setOverrideCursor(Qt::WaitCursor);
      gui->load_mesh(mesh_filename, m_polyhedron, *p_polyhedron_pmaps_bag);
      QApplication::restoreOverrideCursor();

      gui->draw_or_redraw_mesh(
          m_polyhedron,
          p_polyhedron_pmaps_bag,
          viewer_polyhedron,
          false,
          false,
          FEVV::FileUtils::get_file_full_name(mesh_filename));
    }
    FEVV::Block::end("loading-polyhedron");
  }

  ///// Surface_mesh
  FEVV::MeshSurface *m_surface = nullptr;
  if(open_with == OPEN_WITH_SURFACEMESH || open_with == OPEN_WITH_ALL)
  {
    // The following container is already deleted by the viewer destructor
    /// [Snippet Displaying Surface SurfaceMesh]
    FEVV::PMapsContainer *p_surfacemesh_pmaps_bag = new FEVV::PMapsContainer;
    FEVV::Block::begin("loading-surface", "Loading SurfaceMesh mesh.");
    {
      QApplication::setOverrideCursor(Qt::WaitCursor);
      gui->load_mesh(mesh_filename, m_surface, *p_surfacemesh_pmaps_bag);
      QApplication::restoreOverrideCursor();

      gui->draw_or_redraw_mesh(
          m_surface,
          p_surfacemesh_pmaps_bag,
          viewer_surface,
          false,
          false,
          FEVV::FileUtils::get_file_full_name(mesh_filename));
    }
    FEVV::Block::end("loading-surface");
    /// [Snippet Displaying Surface SurfaceMesh]
  }

  ///// Linear_cell_complex
  FEVV::MeshLCC *m_lcc = nullptr;
  if(open_with == OPEN_WITH_LCC || open_with == OPEN_WITH_ALL)
  {
    FEVV::PMapsContainer *p_lcc_pmaps_bag =
        new FEVV::PMapsContainer; // already destroy by the viewer destructor
    FEVV::Block::begin("loading-lcc", "Loading LCC mesh.");
    {
      QApplication::setOverrideCursor(Qt::WaitCursor);
      gui->load_mesh(mesh_filename, m_lcc, *p_lcc_pmaps_bag);
      QApplication::restoreOverrideCursor();

      gui->draw_or_redraw_mesh(
          m_lcc,
          p_lcc_pmaps_bag,
          viewer_lcc,
          false,
          false,
          FEVV::FileUtils::get_file_full_name(mesh_filename));
    }
    FEVV::Block::end("loading-lcc");
  }
#endif // FEVV_USE_CGAL

#ifdef FEVV_USE_OPENMESH
  ///// OpenMesh
  FEVV::MeshOpenMesh *m_openmesh = nullptr;
  if(open_with == OPEN_WITH_OPENMESH || open_with == OPEN_WITH_ALL)
  {
    FEVV::PMapsContainer *p_openmesh_pmaps_bag =
        new FEVV::PMapsContainer; // already destroy by the viewer destructor
    FEVV::Block::begin("loading-openmesh", "Loading OpenMesh mesh.");
    {
      QApplication::setOverrideCursor(Qt::WaitCursor);
      gui->load_mesh(mesh_filename, m_openmesh, *p_openmesh_pmaps_bag);
      QApplication::restoreOverrideCursor();

      gui->draw_or_redraw_mesh(
          m_openmesh,
          p_openmesh_pmaps_bag,
          viewer_openmesh,
          false,
          false,
          FEVV::FileUtils::get_file_full_name(mesh_filename));
    }
    FEVV::Block::end("loading-openmesh");
  }
#endif // FEVV_USE_OPENMESH

#ifdef FEVV_USE_AIF
  ///// AIF
  FEVV::MeshAIF *m_aif = nullptr;
  if(open_with == OPEN_WITH_AIF || open_with == OPEN_WITH_ALL)
  {
#if 1 // AIF nor ready for generic reader
    FEVV::PMapsContainer *p_aif_pmaps_bag =
        new FEVV::PMapsContainer; // already destroy by the viewer destructor
    FEVV::Block::begin("loading-aif", "Loading AIF mesh.");
    {
      QApplication::setOverrideCursor(Qt::WaitCursor);
      gui->load_mesh(mesh_filename, m_aif, *p_aif_pmaps_bag);
      QApplication::restoreOverrideCursor();

      gui->draw_or_redraw_mesh(
          m_aif,
          p_aif_pmaps_bag,
          viewer_aif,
          false,
          false,
          FEVV::FileUtils::get_file_full_name(mesh_filename));
    }
    FEVV::Block::end("loading-aif");
#endif
  }
#endif // FEVV_USE_AIF
  /////// Loading stuff

  // gui->sortModelList();

  int ret = app.exec();

  /////// Cleaning stuff
  if(open_with != OPEN_WITH_NONE)
  {
#ifdef FEVV_USE_CGAL
    // delete viewer_polyhedron;	// already destroy by the adapter
    // destructor delete viewer_surface;		// already destroy by the
    // adapter destructor delete viewer_lcc;		// already destroy by
    // the adapter destructor

    if(!USE_MDI)
      if(adapter_polyhedron != nullptr)
        delete adapter_polyhedron;
    if(!USE_MDI)
      if(adapter_surface != nullptr)
        delete adapter_surface;
    if(!USE_MDI)
      if(adapter_lcc != nullptr)
        delete adapter_lcc;

        // delete m_polyhedron;		// already destroy by the viewer
        // destructor delete m_surface;		// already destroy by the viewer
        // destructor
        // delete m_lcc;			// already destroy by the viewer
        // destructor
#endif
#ifdef FEVV_USE_OPENMESH
    // delete viewer_openmesh;	// already destroy by the adapter destructor
    if(!USE_MDI)
      if(adapter_openmesh != nullptr)
        delete adapter_openmesh;
        // delete m_openmesh;		// already destroy by the viewer
        // destructor
#endif
#ifdef FEVV_USE_AIF
    // delete viewer_aif;		// already destroy by the adapter
    // destructor
    if(!USE_MDI)
      if(adapter_aif != nullptr)
        delete adapter_aif;
        // delete m_aif;			// already destroy by the viewer
        // destructor
#endif
  }

  delete gui;
  /////// Cleaning stuff

  return ret;
}
