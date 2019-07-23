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
// #include "Visualization/SimpleWindow.h"
#include "Visualization/Helpers/QtHelpers.h"

#include <QDebug>
#include <QMenuBar>

#include <QCloseEvent>

#include <boost/assert.hpp>

// Plugins
#include <QApplication> // for qApp
#include <QFileDialog>
#include <QFileInfo>

#if(FEVV_USE_QT5)
#include <QScreen> // for grabWindow
#endif

#include <QPluginLoader>
#include "Visualization/Plugins/PluginInterface.h"
#include "Visualization/Plugins/PluginDialog.hpp"
// Plugins

#ifndef Q_MOC_RUN // MT : very important to avoid the error : ' Parse error at
                  // "BOOST_JOIN" ' -> (qt4 pb with boost)
#include "Visualization/BaseViewerOSG.h"

#include "Visualization/SimpleAdapterVisu.h" // for create an adapter (a child)
#endif


#define USE_MDI true

#define STEP_SPACE 0. // TEMP


#ifndef Q_MOC_RUN // MT : very important to avoid the error : ' Parse error at
                  // "BOOST_JOIN" ' -> (qt4 pb with boost)
#include "FEVV/Filters/Generic/generic_reader.hpp"
#include "FEVV/Filters/Generic/generic_writer.hpp"

#ifdef FEVV_USE_CGAL
#include "FEVV/Filters/CGAL/Point_set/cgal_point_set_reader.hpp"
#include "FEVV/Filters/CGAL/Point_set/cgal_point_set_writer.hpp"
#endif //FEVV_USE_CGAL

#ifdef FEVV_USE_PCL
#include "FEVV/Filters/PCL/pcl_point_cloud_reader.hpp"
#include "FEVV/Filters/PCL/pcl_point_cloud_writer.hpp"
#endif //FEVV_USE_PCL
#endif //Q_MOC_RUN


inline FEVV::SimpleWindow::SimpleWindow(QWidget *_parent,
                                        Qt::WindowFlags _flags)
    : BaseWindowQt(_parent, _flags)
{
#ifdef DEBUG_VISU2
  std::cout << "*** this=" << this << "    entering " << __func__ << std::endl;
#endif

  ui.setupUi(this);

#ifdef DEBUG_VISU2
  std::cout << "*** this=" << this << "    leaving " << __func__ << std::endl;
#endif
}

inline FEVV::SimpleWindow::~SimpleWindow()
{
#ifdef DEBUG_VISU2
  std::cout << "*** this=" << this << "    entering " << __func__ << std::endl;
#endif

  if(mdiArea != nullptr)
  {
    delete mdiArea;
  }

#ifdef DEBUG_VISU2
  std::cout << "*** this=" << this << "    leaving " << __func__ << std::endl;
#endif
}

inline void
FEVV::SimpleWindow::attach(Adapter *_adapter)
{
  AdapterQt *adapter = dynamic_cast< AdapterQt * >(_adapter);
  if(adapter)
  {
    attach(adapter);
  }
  else
  {
    Assert::check(false,
                  "is not implemented. See attach(AdapterQt) instead.",
                  "SimpleWindow::attach(Adapter)");
  }
}

inline void
FEVV::SimpleWindow::attach(AdapterQt *_adapter, const bool _useMdiWindows)
{
  if(!Assert::check(_adapter != nullptr,
                    "The given adapter is null.",
                    "SimpleWindow::attach(AdapterQt,bool)"))
  {
    return;
  }

  adapters.push_back(_adapter);
  _adapter->getViewer()->attach(this);

  if(_useMdiWindows)
  {
    if(mdiArea == nullptr) // Only once.
    {
      mdiArea = new QMdiArea(this);
      // mdiArea->setHorizontalScrollBarPolicy( Qt::ScrollBarAsNeeded ); //
      // Qt::ScrollBarAlwaysOff mdiArea->setVerticalScrollBarPolicy(
      // Qt::ScrollBarAsNeeded ); // Qt::ScrollBarAlwaysOff

      ui.gridLayout->addWidget(mdiArea, 0, 0);

      connect(
          ui.actionTile, SIGNAL(triggered()), mdiArea, SLOT(tileSubWindows()));
      connect(ui.actionCascade,
              SIGNAL(triggered()),
              mdiArea,
              SLOT(cascadeSubWindows()));

      // connect(mdiArea, SIGNAL(subWindowActivated(QMdiSubWindow *)), this,
      // SLOT(updateMenus())); // TODO
    }

    _adapter->setMinimumSize(300, 200); // default value is 300 x 200 pixels
    mdiArea->addSubWindow(_adapter);
  }
  else
  {
    unsigned int nbAdapters = adapters.size();
    unsigned int nbLines = (unsigned int)round(std::sqrt(nbAdapters));
    unsigned int nbColumns = (unsigned int)ceil((double)nbAdapters / nbLines);

    for(unsigned int iLines = 0; iLines < nbLines; ++iLines)
    {
      for(unsigned int iColumns = 0; iColumns < nbColumns; ++iColumns)
      {
        unsigned int index = iLines * nbColumns + iColumns;
        if(index >= nbAdapters)
        {
          break;
        }
        ui.gridLayout->addWidget(
            static_cast< AdapterQt * >(adapters[index]), iLines, iColumns);

        // if( index % 2 == 0 )
        // {
        // ui.gridLayout->addWidget( static_cast<AdapterQt*>(adapters[index]),
        // index, 0 );//iLines+1, iColumns+1 );
        // }
        // else
        // {
        // ui.gridLayout->addWidget( static_cast<AdapterQt*>(adapters[index]),
        // 0, index );//iLines+1, iColumns+1 );
        // }
      }
    }
  }

  update();
}

inline void
FEVV::SimpleWindow::attachPlugin(Plugin *_plugin)
{
  _plugin->addParameters(this);
}

inline void
FEVV::SimpleWindow::init()
{
  init(false);
}

inline void
FEVV::SimpleWindow::init(const bool _test, const int _width, const int _height)
{
  if(!Assert::check(
         !bIsInit, "is already init. Leaving...", "SimpleWindow::init"))
  {
    return;
  }

  setlocale(LC_ALL, "C"); // very important for Linux

  // --- JUST HERE FOR AUTOMATIC TEST ---
  if(_test)
  {
    connect(&timerQuit, SIGNAL(timeout()), this, SLOT(close()));
    timerQuit.start(8000);
  }
  // --- JUST HERE FOR AUTOMATIC TEST ---

#ifdef DEBUG_VISU
  Helpers::changeBackgroundColor(this, Color::Green());
  Helpers::changeBackgroundColor(ui.statusbar, Color::Pin());
  Helpers::changeBackgroundColor(ui.dockWidget, Color::Blue());
#endif

  setMinimumSize(_width, _height);
  resize(_width, _height);

  std::srand(std::time(0)); // @todo only for debug

  QObject::connect(
      ui.listModels,
      SIGNAL(currentItemChanged(QListWidgetItem *, QListWidgetItem *)),
      this,
      SLOT(onItemChanged(QListWidgetItem *, QListWidgetItem *)));

  QMenuBar *menuBar = this->menuBar();

  // QMenu* menu = menuBar->addMenu( "Debug" );
  // menu->addAction( "Add Random Sphere", this, SLOT( onAddBall() ) );

  // Plugins
  menuPlugins = ui.menuPlugins; // menuBar->addMenu("Plugins (old)");
  menuPlugins->addAction("Plugins information", this, SLOT(aboutPlugins()));

#if(FEVV_USE_QT5)
  ui.menuHelp->addSeparator();
  ui.menuHelp->addAction("Screenshot", this, SLOT(onGrab()));
#endif

  QObject::connect(
      ui.applyButton, SIGNAL(clicked(bool)), this, SLOT(onApplyButton()));
  ui.applyButton->setEnabled(false);

  ui.listParams->setVisible(false);

  // time
  // TODO LATER - SPACE/TIME - FINDME
  /*ui.actionDyn_First->setEnabled(false);
  ui.actionDyn_Previous->setEnabled(false);
  ui.actionDyn_Next->setEnabled(false);
  ui.actionDyn_Last->setEnabled(false);

  ui.actionChange_viewer_mode->setEnabled(false);*/
  // time

  // create default empty viewer
  //createNewViewer();

  bIsInit = true;
}

inline void
FEVV::SimpleWindow::setParam(std::string _name,
                             int *_value,
                             std::string _pluginName,
                             Plugin *_plugin)
{
#if 0
    // @todo missing delete of next line allocation
    QLabel* testLabel = new QLabel( ui.listParams );
    testLabel->setObjectName(QString::fromUtf8("label"));
    testLabel->setText(QString::fromUtf8(_name.c_str()));

    // @todo missing delete of next line allocation
    SimpleLineEdit* testValue = new SimpleLineEdit( ui.listParams );

    testValue->setObjectName(QString::fromUtf8("value"));
    testValue->setText(QString::fromUtf8(std::to_string(*_value).c_str()));
    testValue->setParams(_value, _pluginName, _plugin);

    QObject::connect(testValue, SIGNAL(textChanged(QString)), testValue, SLOT(modificationSlot()));
    QObject::connect(testValue, SIGNAL(modificationSignal(std::string, BasePlugin*)), this, SLOT(onModificationParam(std::string, BasePlugin*)));

	QObject* p = dynamic_cast< QObject* >(_plugin);
	if (p) QObject::connect(p, SIGNAL(resetSignal()), testValue, SLOT(resetSlot()));

    ui.formLayout->addRow( testLabel, testValue );
#endif
}

inline void
FEVV::SimpleWindow::setParam(std::string _name,
                             double *_value,
                             std::string _pluginName,
                             Plugin *_plugin)
{
#if 0
    QLabel* testLabel = new QLabel( ui.listParams );
    testLabel->setObjectName(QString::fromUtf8("label"));
    testLabel->setText(QString::fromUtf8(_name.c_str()));

    SimpleLineEdit* testValue = new SimpleLineEdit( ui.listParams );
    testValue->setObjectName(QString::fromUtf8("value"));
    testValue->setText(QString::fromUtf8(std::to_string(*_value).c_str()));
    testValue->setParams(_value, _pluginName, _plugin);

    QObject::connect(testValue, SIGNAL(textChanged(QString)), testValue, SLOT(modificationSlot()));
    QObject::connect(testValue, SIGNAL(modificationSignal(std::string, BasePlugin*)), this, SLOT(onModificationParam(std::string, BasePlugin*)));

	QObject* p = dynamic_cast< QObject* >(_plugin);
	if (p) QObject::connect(p, SIGNAL(resetSignal()), testValue, SLOT(resetSlot()));

    ui.formLayout->addRow( testLabel, testValue );
#endif
}

inline void
FEVV::SimpleWindow::setParam(std::string _name,
                             float *_value,
                             std::string _pluginName,
                             Plugin *_plugin)
{
#if 0
    QLabel* testLabel = new QLabel( ui.listParams );
    testLabel->setObjectName(QString::fromUtf8("label"));
    testLabel->setText(QString::fromUtf8(_name.c_str()));

    SimpleLineEdit* testValue = new SimpleLineEdit( ui.listParams );
    testValue->setObjectName(QString::fromUtf8("value"));
    testValue->setText(QString::fromUtf8(std::to_string(*_value).c_str()));
    testValue->setParams(_value, _pluginName, _plugin);

    QObject::connect(testValue, SIGNAL(textChanged(QString)), testValue, SLOT(modificationSlot()));
    QObject::connect(testValue, SIGNAL(modificationSignal(std::string, BasePlugin*)), this, SLOT(onModificationParam(std::string, BasePlugin*)));

	QObject* p = dynamic_cast< QObject* >(_plugin);
	if (p) QObject::connect(p, SIGNAL(resetSignal()), testValue, SLOT(resetSlot()));

    ui.formLayout->addRow( testLabel, testValue );
#endif
}

inline void
FEVV::SimpleWindow::setParam(std::string _name,
                             bool *_value,
                             std::string _pluginName,
                             Plugin *_plugin)
{
#if 0
    QLabel* testLabel = new QLabel( ui.listParams );
    testLabel->setObjectName(QString::fromUtf8("label"));
    testLabel->setText(QString::fromUtf8(_name.c_str()));

    SimpleCheckBox* testValue = new SimpleCheckBox( ui.listParams );
    testValue->setObjectName(QString::fromUtf8("value"));
    Qt::CheckState state = (*_value) ? Qt::CheckState::Checked : Qt::CheckState::Unchecked;
    testValue->setCheckState(state);
    testValue->setParams(_value, _pluginName, _plugin);

    QObject::connect(testValue, SIGNAL(stateChanged(int)), testValue, SLOT(modificationSlot()));
    QObject::connect(testValue, SIGNAL(modificationSignal(std::string, BasePlugin*)), this, SLOT(onModificationParam(std::string, BasePlugin*)));

	QObject* p = dynamic_cast< QObject* >(_plugin);
	if (p) QObject::connect(p, SIGNAL(resetSignal()), testValue, SLOT(resetSlot()));

    ui.formLayout->addRow( testLabel, testValue );
#endif
}

inline void
FEVV::SimpleWindow::setParam(std::string _name,
                             std::string *_value,
                             std::string _pluginName,
                             Plugin *_plugin)
{
#if 0
    QLabel* testLabel = new QLabel( ui.listParams );
    testLabel->setObjectName(QString::fromUtf8("label"));
    testLabel->setText(QString::fromUtf8(_name.c_str()));

    SimpleLineEdit* testValue = new SimpleLineEdit( ui.listParams );
    testValue->setObjectName(QString::fromUtf8("value"));
    testValue->setText(QString::fromUtf8(_value->c_str()));
    testValue->setParams(_value, _pluginName, _plugin);

    QObject::connect(testValue, SIGNAL(textChanged(QString)), testValue, SLOT(modificationSlot()));
    QObject::connect(testValue, SIGNAL(modificationSignal(std::string, BasePlugin*)), this, SLOT(onModificationParam(std::string, BasePlugin*)));

	QObject* p = dynamic_cast< QObject* >(_plugin);
	if (p) QObject::connect(p, SIGNAL(resetSignal()), testValue, SLOT(resetSlot()));

    ui.formLayout->addRow( testLabel, testValue );
#endif
}

inline void
FEVV::SimpleWindow::onModificationParam(std::string _pluginName,
                                        Plugin *_plugin)
{
  // std::cout << "onModificationParam called" << std::endl;
  if(_pluginName == "")
  {
    return;
  }
  if(stackPlugins.find(_pluginName) == stackPlugins.end())
  {
    stackPlugins.insert({_pluginName, _plugin});
    ui.applyButton->setEnabled(true);
  }
  // else
  // {
  //     std::cout << "Already inside stackPlugins" << std::endl;
  // }
}

inline void
FEVV::SimpleWindow::onApplyButton()
{
  // if no viewer is opened, open one
  if(!isValid())
    auto viewer = createNewViewer();

  // apply plugin in current viewer
  for(auto funct : stackPlugins)
  {
    for(auto w: adapters)
    {
      if(w->isSelected())
        w->apply(funct.second);
    }
  }

  stackPlugins.clear();
  ui.applyButton->setEnabled(false);
}

// Plugins
inline void
FEVV::SimpleWindow::loadQtPlugins()
{
  foreach(QObject *plugin, QPluginLoader::staticInstances())
  {
    populateMenus(plugin);
  }

  // std::cout << "-> loadPlugins from " <<
  // /*qApp->*/QApplication::applicationDirPath().toStdString() << "\n";
  pluginsDir = QDir(/*qApp->*/ QApplication::applicationDirPath());

#if defined(Q_OS_WIN)
  /*if (pluginsDir.dirName().toLower() == "debug" ||
  pluginsDir.dirName().toLower() == "release") pluginsDir.cdUp();*/
#elif defined(Q_OS_MAC)
  if(pluginsDir.dirName() == "MacOS")
  {
    pluginsDir.cdUp();
    pluginsDir.cdUp();
    pluginsDir.cdUp();
  }
#endif
  // pluginsDir.cd("plugins");

  foreach(QString fileName, pluginsDir.entryList(QDir::Files))
  {
    if(/*fileName.contains("Filter") && */ QLibrary::isLibrary(fileName))
    {
      QPluginLoader loader(pluginsDir.absoluteFilePath(fileName));
      QObject *plugin = loader.instance();
      if(plugin)
      {
        menuPlugins->addSeparator();
        populateMenus(plugin);
        pluginFileNames += fileName;

        // call plugin's init() methods
        Generic_PluginInterface *mepp_qt_plugin =
            qobject_cast< Generic_PluginInterface * >(plugin);
        mepp_qt_plugin->init(this);

        FEVV::BasePlugin *mepp_plugin =
            dynamic_cast< FEVV::BasePlugin * >(plugin);
        mepp_plugin->init();
        this->attachPlugin(mepp_plugin);
        // call plugin's init() methods
      }
    }
  }

  menuPlugins->setEnabled(!menuPlugins->actions().isEmpty());
}

inline void
FEVV::SimpleWindow::populateMenus(QObject *plugin)
{
  Generic_PluginInterface *iGeneric_Filter =
      qobject_cast< Generic_PluginInterface * >(plugin);
  if(iGeneric_Filter)
  {
    addToMenu(plugin,
              iGeneric_Filter->Generic_plugins(),
              menuPlugins,
              SLOT(applyPlugin()));
  }
}

inline void
FEVV::SimpleWindow::addToMenu(QObject *plugin,
                              const QStringList &texts,
                              QMenu *menu,
                              const char *member,
                              QActionGroup *actionGroup)
{
  foreach(QString text, texts)
  {
    QAction *action = new QAction(text, plugin);
    connect(action, SIGNAL(triggered()), this, member);
    menu->addAction(action);

    if(actionGroup)
    {
      action->setCheckable(true);
      actionGroup->addAction(action);
    }
  }
}

inline void
FEVV::SimpleWindow::applyPlugin()
{
  QAction *action = qobject_cast< QAction * >(sender());

  Generic_PluginInterface *iGeneric_Filter =
      qobject_cast< Generic_PluginInterface * >(action->parent());
  if(iGeneric_Filter)
  {
    iGeneric_Filter->Generic_plugin(action->text());
  }
}

inline void
FEVV::SimpleWindow::aboutPlugins()
{
  PluginDialog dialog(pluginsDir.path(), pluginFileNames, this);
  dialog.exec();
}
// Plugins

inline void
FEVV::SimpleWindow::sortModelList()
{
  ui.listModels->sortItems(Qt::AscendingOrder); // Qt::DescendingOrder
}

inline void
FEVV::SimpleWindow::notify()
{
  if(!Assert::check(isValid(),
                    "is not valid (see init() or attach()). Leaving...",
                    "SimpleWindow::notify"))
  {
    return;
  }

  update();
}

inline void
FEVV::SimpleWindow::update(bool pick)
{
  updateModelList(pick);
}

inline void
FEVV::SimpleWindow::enableSpaceTimeMenus()
{
  // TODO LATER - SPACE/TIME - FINDME
  /*ui.actionChange_viewer_mode->setEnabled(true);

  ui.actionChange_viewer_mode->setText(QObject::tr("Change viewer mode (-> to
  TIME mode)"));*/
}

inline void
FEVV::SimpleWindow::updateModelList(bool pick)
{
  pick = true;
  //ELO-note: forcing the pick flag to true make mesh selection work as
  //          expected in SPACE/TIME mode. Don't know why.
  using Viewer = BaseViewerOSG;

  ui.listModels->clear(); //<! All items will be permanently deleted.
  if(adapters.size() > 0)
  {
    for(unsigned int aa = 0; aa < adapters.size(); ++aa)
    {
      QListWidgetItem *pickitem = nullptr;

      Viewer *viewer = static_cast< Viewer * >(adapters[aa]->getViewer());
      Viewer::DataModelVector *result = viewer->getDataModel();

      if(result->size() > 0)
      {
        for(unsigned int rr = 0; rr < result->size(); ++rr)
        {
          if((*result)[rr].type == Helpers::DataType::MODEL)
          {
            QListWidgetItem *item = new QListWidgetItem(
                QString::fromUtf8((*result)[rr].name.c_str()));
            QVariant variant;
            variant.setValue((*result)[rr]);
            item->setData(Qt::UserRole, variant);

            ui.listModels->insertItem(0, item); // addItem(item);

            // --------------------------------------------------------------------------------
            if(pick)
            {
              // if ( (*result)[rr].node != nullptr ) // not necessary ?
              {
                // std::cout << "UPDATE : " << (*result)[rr].node->getName() <<
                // " - " << (*result)[rr].node << std::endl;

                if((*result)[rr].viewer->isNodeSelected(
                       (*result)[rr]
                           .node)) // [PICK] : IMPORTANT test in PickHandler.h
                                   // --> see keyword [PICK] - (if we re-select
                                   // same node, HERE the function doesn't
                                   // return TRUE ! -> why ?)
                {
                  pickitem = item;
                  // std::cout << "pickitem" << std::endl;
                }
              }
            }
            // --------------------------------------------------------------------------------
          }
        }

        // ----------------------------------------------------------------------------------------
        if(pick)
        {
          if(pickitem != nullptr)
          {
            // std::cout << "setCurrentItem pickitem" << std::endl;
            ui.listModels->setCurrentItem(pickitem);
          }
        }
        // ----------------------------------------------------------------------------------------
        else
          ui.listModels->setCurrentRow(0);

        // QMdiSubWindow FOCUS
        /*if (mdiArea) // COMMENTED 01/12/17
        {
                  BaseAdapterVisuQt *bavQt = dynamic_cast< BaseAdapterVisuQt*
        >(adapters[aa]); QList<QMdiSubWindow *> listMdiSubW =
        mdiArea->subWindowList();

                  for (int MdiSubW = 0; MdiSubW < listMdiSubW.size(); ++MdiSubW)
                  {
                          BaseAdapterVisuQt *bavQtMdiSubW = dynamic_cast<
        BaseAdapterVisuQt* >(listMdiSubW.at(MdiSubW)->widget()); if
        (bavQtMdiSubW->getViewer() == bavQt->getViewer())
                          {
                                mdiArea->setActiveSubWindow(listMdiSubW.at(MdiSubW));
                                //std::cout << "updateModelList(): " <<
        listMdiSubW.at(MdiSubW)->windowTitle().toStdString() << std::endl;
                                break;
                          }
                  }
        }*/
        // QMdiSubWindow FOCUS
      }
    }
  }
}

template< typename MeshT >
void
FEVV::SimpleWindow::draw_or_redraw_mesh(
    /*const */ MeshT *mesh,
    /*const */ FEVV::PMapsContainer *pmaps_bag,
    FEVV::SimpleViewer *viewer,
    bool _redraw,
    bool _recomputeNT_if_redraw,
    std::string _mesh_filename,
    bool _recreateOSGobj_if_redraw,
    float _step)
{
  // QApplication::setOverrideCursor(Qt::BusyCursor);
  viewer->draw_or_redraw_mesh(mesh,
                              pmaps_bag,
                              _redraw,
                              _recomputeNT_if_redraw,
                              _mesh_filename,
                              _recreateOSGobj_if_redraw,
                              _step);
  // QApplication::restoreOverrideCursor();
}

inline void
FEVV::SimpleWindow::on_actionNew_triggered()
{
  QMessageBox::information(this, "", QObject::tr("New - Not yet implemented."));
}


template< typename HalfedgeGraph >
inline void
FEVV::SimpleWindow::on_actionOpen_SPACE_TIME(FEVV::SimpleViewer *viewer)
{
  std::string ds_name =
      FEVV::getDatastructureName(static_cast< HalfedgeGraph * >(nullptr));

  QString allExtensions;

  if(ds_name == "CGALPOINTSET")
  {
    allExtensions = "XYZ/OFF/PLY files (*.xyz *.off *.ply);;"
                    "XYZ files (*.xyz);;"
                    "OFF files (*.off);;"
                    "PLY files (*.ply)";
  }
  else if(ds_name == "PCLPOINTCLOUD")
  {
    allExtensions = "XYZ/PCD/PLY files (*.xyz *.pcd *.ply);;"
                    "XYZ files (*.xyz);;"
                    "PCD files (*.pcd);;"
                    "PLY files (*.ply)";
  }
  else
  {
    QString defaultExtensions = "OBJ/OFF files (*.obj *.off);;"
                                "OBJ files (*.obj);;"
                                "OFF files (*.off);;"
                                "COFF files (*.coff);;"
                                "PLY files (*.ply);;"
                                "MSH files (*.msh)";

    QString vtkExtensions = "VTK Files (*.vtk);;"
                            "VTP files (*.vtp);;"
                            "VTU files (*.vtu)";

    allExtensions = defaultExtensions;
#ifdef FEVV_USE_VTK
    allExtensions += ";;" + vtkExtensions;
#endif
#ifdef FEVV_USE_FBX
    allExtensions += ";;FBX files (*.fbx)";
#endif
  }

  QString suffix;
  QFileDialog::Option options = (QFileDialog::Option)0;

#if defined(__linux__) || defined(__APPLE__)
  options = QFileDialog::DontUseNativeDialog; // PB under LINUX !?
#endif

  QStringList files_qt =
      QFileDialog::getOpenFileNames(this,
                                    "Open (SPACE/TIME)",
                                    /*openLocation*/ QDir::currentPath(),
                                    allExtensions,
                                    &suffix,
                                    options);

  // convert QStringList to standard type
  std::vector< std::string > files;
  for(auto qstr: files_qt)
    files.push_back(qstr.toStdString());

  // open files
  open_SPACE_TIME< HalfedgeGraph >(viewer, files);
}


template< typename HalfedgeGraph >
inline void
FEVV::SimpleWindow::open_SPACE_TIME(FEVV::SimpleViewer *viewer,
                                    const std::vector< std::string >& files)
{
  // load and draw meshes
  int m = 0;
  for(auto filename: files)
  {
    HalfedgeGraph *mesh = new HalfedgeGraph;
      // destroyed by the viewer destructor
    FEVV::PMapsContainer *p_pmaps_bag = new FEVV::PMapsContainer;
      // destroyed by the viewer destructor

    QApplication::setOverrideCursor(Qt::WaitCursor);
    statusBar()->showMessage(QObject::tr("Open mesh file...") /*, 2000*/);

    try
    {
      // read mesh from file
      FEVV::Filters::read_mesh(
          filename, *mesh, *p_pmaps_bag, open_only_pts_mode);

      // open a new viewer if needed
      // and only if the mesh reading was successfull (aka no exception)
      // to avoid empty viewer
      if(viewer == nullptr || ctrl_pressed) // NOTE : ctrl_pressed not documented (here only for internal usage...)
        viewer = createNewViewer();

      // draw mesh
      draw_or_redraw_mesh(
          mesh,
          p_pmaps_bag,
          viewer,
          false,
          false,
          FEVV::FileUtils::get_file_full_name(filename),
          true,
          m * STEP_SPACE);

      ++m;
    }
    catch(const std::exception& e)
    {
      std::cout << e.what() << '\n';
      QMessageBox::warning(
          0, "", QObject::tr(e.what()));

      // reading the mesh failed, so release memory now
      delete mesh;
      delete p_pmaps_bag;
    }

    updateActiveChildTitle();

    statusBar()->showMessage(QObject::tr("") /*, 2000*/);
    QApplication::restoreOverrideCursor();
  }
}


inline void
FEVV::SimpleWindow::on_actionOpen_triggered()
{
  // capture keyboard state
  bool shift_pressed =
      QApplication::keyboardModifiers().testFlag(Qt::ShiftModifier);
  bool alt_pressed =
      QApplication::keyboardModifiers().testFlag(Qt::AltModifier);
  ctrl_pressed =
      QApplication::keyboardModifiers().testFlag(Qt::ControlModifier);

  // for 'only_pts' mode
  if(alt_pressed)
    open_only_pts_mode = true;
  else
    open_only_pts_mode = false;

  // logic: default-mode is "add in SPACE/TIME" mode in current viewer ;
  //        if no viewer is opened, then open one in SPACE mode

  // use existing viewer unless shift is pressed
  // note: the call to activeMdiChild() to retrieve the current viewer 
  //       must come before the call to chooseDatastructureMsgBox()
  //       because the latter resets the active window
  SimpleViewer *viewer = nullptr;
  if( activeMdiChild() && (!shift_pressed) && (!ctrl_pressed) ) // NOTE : ctrl_pressed not documented (here only for internal usage...)
  {
    // open mesh in current viewer
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());
    Adapter::Viewer *base_viewer = bavQt->getViewer();
    viewer = dynamic_cast< SimpleViewer * >(base_viewer);
    assert(viewer != nullptr);
  }

  // ask the user for the datastructure
  std::string mesh_type = chooseDatastructureMsgBox();
  if(mesh_type == "NONE")
    return; // cancel pressed, aborting

  // open mesh(es)
#ifdef FEVV_USE_CGAL
  if(mesh_type == "POLYHEDRON")
    on_actionOpen_SPACE_TIME< FEVV::MeshPolyhedron >(viewer);
  else if(mesh_type == "SURFACEMESH")
    on_actionOpen_SPACE_TIME< FEVV::MeshSurface >(viewer);
  else if(mesh_type == "LCC")
    on_actionOpen_SPACE_TIME< FEVV::MeshLCC >(viewer);
  else if(mesh_type == "CGALPOINTSET")
    on_actionOpen_SPACE_TIME< FEVV::CGALPointSet >(viewer);
#endif

#ifdef FEVV_USE_OPENMESH
  if(mesh_type == "OPENMESH")
    on_actionOpen_SPACE_TIME< FEVV::MeshOpenMesh >(viewer);
#endif

#ifdef FEVV_USE_AIF
  if(mesh_type == "AIF")
    on_actionOpen_SPACE_TIME< FEVV::MeshAIF >(viewer);
#endif

#ifdef FEVV_USE_PCL
  if(mesh_type == "PCLPOINTCLOUD")
    on_actionOpen_SPACE_TIME< FEVV::PCLPointCloud >(viewer);
#endif
}


inline
void
FEVV::SimpleWindow::writeHG(FEVV::SimpleViewer *viewer)
{
  if(!viewer)
    return;

  QString defaultExtensions = "OBJ files (*.obj);;"
                              "OFF files (*.off);;"
                              "COFF files (*.coff);;"
                              "PLY files (*.ply);;"
                              "MSH files (*.msh)";

  QString vtkExtensions = "VTK files (*.vtk);;"
                          "VTP files (*.vtp);;"
                          "VTU files (*.vtu)";

  QString cgalpointsetExtensions = "XYZ files (*.xyz);;"
                                   "OFF files (*.off);;"
                                   "PLY files (*.ply)";

  QString pclpointcloudExtensions = "PCD files (*.pcd);;"
                                    "PLY files (*.ply)";

  FEVV::MixedMeshesVector meshes = viewer->getSelectedMeshes();
  std::vector< std::string > meshes_names = viewer->getSelectedMeshesNames();
  std::vector< FEVV::PMapsContainer * > properties_maps =
      viewer->getSelected_properties_maps();

  for(unsigned i = 0; i < meshes.size(); i++)
  {
    QString allExtensions;

    if(meshes[i].second == "CGALPOINTSET")
    {
      allExtensions = cgalpointsetExtensions;
    }
    else if(meshes[i].second == "PCLPOINTCLOUD")
    {
      allExtensions = pclpointcloudExtensions;
    }
    else
    {
      allExtensions = defaultExtensions;
#ifdef FEVV_USE_VTK
      allExtensions += ";;" + vtkExtensions;
#endif
    }
    

    QString suffix;
    QFileDialog::Option options = (QFileDialog::Option)0;

#ifdef __linux__
    options = QFileDialog::DontUseNativeDialog; // PB under LINUX !?
#endif

    QString fileName = QFileDialog::getSaveFileName(
        0,
        "Save As",
        /*saveLocation*/ /*QDir::currentPath()*/
        QFileInfo(QString::fromStdString(meshes_names[i]))
            .baseName(),
        allExtensions,
        &suffix,
        options);

#ifdef __linux__
    if(suffix.indexOf(".obj") >= 0)
      fileName += ".obj";
    else if(suffix.indexOf(".off") >= 0)
      fileName += ".off";
    else if(suffix.indexOf(".coff") >= 0)
      fileName += ".coff";
    else if(suffix.indexOf(".ply") >= 0)
      fileName += ".ply";
    else if(suffix.indexOf(".msh") >= 0)
      fileName += ".msh";
    else if(suffix.indexOf(".vtk") >= 0)
      fileName += ".vtk";
    else if(suffix.indexOf(".vtp") >= 0)
      fileName += ".vtp";
    else if(suffix.indexOf(".vtu") >= 0)
      fileName += ".vtu";
    else if(suffix.indexOf(".xyz") >= 0)
      fileName += ".xyz";
    else if(suffix.indexOf(".pcd") >= 0)
      fileName += ".pcd";
#endif

    if(!fileName.isEmpty())
    {
#ifdef FEVV_USE_CGAL
      if(meshes[i].second == "POLYHEDRON")
      {
        auto mesh_ptr = static_cast< FEVV::MeshPolyhedron* >(meshes[i].first);
        FEVV::Filters::write_mesh(fileName.toStdString(),
                                  *mesh_ptr,
                                  *(properties_maps[i]));
      }
      if(meshes[i].second == "SURFACEMESH")
      {
        auto mesh_ptr = static_cast< FEVV::MeshSurface* >(meshes[i].first);
        FEVV::Filters::write_mesh(fileName.toStdString(),
                                  *mesh_ptr,
                                  *(properties_maps[i]));
      }
      if(meshes[i].second == "LCC")
      {
        auto mesh_ptr = static_cast< FEVV::MeshLCC* >(meshes[i].first);
        FEVV::Filters::write_mesh(fileName.toStdString(),
                                  *mesh_ptr,
                                  *(properties_maps[i]));
      }
      if(meshes[i].second == "CGALPOINTSET")
      {
        auto mesh_ptr = static_cast< FEVV::CGALPointSet* >(meshes[i].first);
        FEVV::Filters::write_mesh(fileName.toStdString(),
                                  *mesh_ptr,
                                  *(properties_maps[i]));
      }
#endif //FEVV_USE_CGAL

#ifdef FEVV_USE_OPENMESH
      if(meshes[i].second == "OPENMESH")
      {
        auto mesh_ptr = static_cast< FEVV::MeshOpenMesh* >(meshes[i].first);
        FEVV::Filters::write_mesh(fileName.toStdString(),
                                  *mesh_ptr,
                                  *(properties_maps[i]));
      }
#endif //FEVV_USE_OPENMESH

#ifdef FEVV_USE_AIF
      if(meshes[i].second == "AIF")
      {
        auto mesh_ptr = static_cast< FEVV::MeshAIF* >(meshes[i].first);
        FEVV::Filters::write_mesh(fileName.toStdString(),
                                  *mesh_ptr,
                                  *(properties_maps[i]));
      }
#endif //FEVV_USE_AIF

#ifdef FEVV_USE_PCL
      if(meshes[i].second == "PCLPOINTCLOUD")
      {
        auto mesh_ptr = static_cast< FEVV::PCLPointCloud* >(meshes[i].first);
        FEVV::Filters::write_mesh(fileName.toStdString(),
                                  *mesh_ptr,
                                  *(properties_maps[i]));
      }
#endif //FEVV_USE_PCL
    }
  }
}


inline void
FEVV::SimpleWindow::on_actionSaveAs_triggered()
{
  for(unsigned i = 0; i < adapters.size(); i++)
  {
    if(adapters[i]->isSelected())
    {
      BaseAdapterVisu *bav = adapters[i];
      auto viewer = dynamic_cast< SimpleViewer * >(bav->getViewer());
      writeHG(viewer);
    }
  }
}


inline void
FEVV::SimpleWindow::on_actionClose_triggered()
{
  if(mdiArea)
    mdiArea->closeActiveSubWindow();
}

inline void
FEVV::SimpleWindow::on_actionQuit_triggered()
{
  close();
}

inline void
FEVV::SimpleWindow::closeEvent(QCloseEvent *event)
{
  if(mdiArea)
  {
    mdiArea->closeAllSubWindows();

    if(activeMdiChild())
    {
      event->ignore();
    }
    else
    {
      event->accept();
    }
  }
  else
  {
    event->accept();
  }
}

inline void
FEVV::SimpleWindow::on_actionClose_window_triggered()
{
  on_actionClose_triggered();
}

inline void
FEVV::SimpleWindow::on_actionClose_all_triggered()
{
  if(mdiArea)
    mdiArea->closeAllSubWindows();
}

inline void
FEVV::SimpleWindow::updateActiveChildTitle()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

    RenderMode Render = bavQt->getViewer()->m_RenderMode;
    unsigned int Mode = (unsigned int)bavQt->getViewer()->m_space_time +
                        (unsigned int)bavQt->getViewer()->m_time;

    QString sRender, sMode;

    if(Render == RenderMode::RENDER_LEGACY)
      sRender = "LEGACY";
    else if(Render == RenderMode::RENDER_SHADERS_DIRECT_LIGHTING)
      sRender = "SHADERS (DIRECT LIGHTING)";
    else if(Render == RenderMode::RENDER_SHADERS_INDIRECT_LIGHTING)
      sRender = "SHADERS (INDIRECT LIGHTING)";

    if(Mode == 0)
      sMode = "NORMAL";
    else if(Mode == 1)
      sMode = "SPACE";
    else if(Mode == 2)
      sMode = "TIME";

    bavQt->setWindowTitle(
        bavQt->windowTitle().left(bavQt->windowTitle().indexOf('>') + 1) +
        QObject::tr(" - Mode : %1 - Render : %2").arg(sMode).arg(sRender));
  }
}

inline void
FEVV::SimpleWindow::on_actionChange_MDI_view_mode_triggered()
{
  if(mdiArea)
  {
    if(mdiArea->viewMode() == QMdiArea::SubWindowView)
    {
      mdiArea->setViewMode(QMdiArea::TabbedView);
      ui.actionChange_MDI_view_mode->setText(
          QObject::tr("Change MDI view mode (-> to subwindow view)"));
    }
    else
    {
      mdiArea->setViewMode(QMdiArea::SubWindowView);
      ui.actionChange_MDI_view_mode->setText(
          QObject::tr("Change MDI view mode (-> to tabbed view)"));
    }
  }
}

inline void
FEVV::SimpleWindow::on_actionChange_viewer_mode_triggered()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

    if(bavQt->getViewer()->m_space_time)
    {
      // TODO LATER - SPACE/TIME - FINDME
      /*if (!(bavQt->getViewer()->m_time))
      {
          ui.actionDyn_First->setEnabled(true);
          ui.actionDyn_Previous->setEnabled(true);
          ui.actionDyn_Next->setEnabled(true);
          ui.actionDyn_Last->setEnabled(true);

          ui.actionChange_viewer_mode->setText(QObject::tr("Change viewer mode
      (-> to SPACE mode)"));
      }
      else
      {
          ui.actionDyn_First->setEnabled(false);
          ui.actionDyn_Previous->setEnabled(false);
          ui.actionDyn_Next->setEnabled(false);
          ui.actionDyn_Last->setEnabled(false);

          ui.actionChange_viewer_mode->setText(QObject::tr("Change viewer mode
      (-> to TIME mode)"));
      }*/

      bavQt->getViewer()->m_time = !(bavQt->getViewer()->m_time);
      // bavQt->setWindowTitle(QObject::tr("MODE:
      // %1").arg(bavQt->getViewer()->m_time));

      // bavQt->getViewer()->m_space_time_changeColorMode = false; // useful or
      // not ?
      pre_actionHG(bavQt->getViewer(), 'M');
      // bavQt->getViewer()->m_space_time_changeColorMode = true;  // useful or
      // not ?

      updateActiveChildTitle();
    }
  }
}

inline void
FEVV::SimpleWindow::activate_time_mode()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

    if(bavQt->getViewer()->m_space_time)
    {
      bavQt->getViewer()->m_time = true;

      pre_actionHG(bavQt->getViewer(), 'M');

      updateActiveChildTitle();
    }
  }
}

inline void
FEVV::SimpleWindow::activate_space_mode()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

    if(bavQt->getViewer()->m_space_time)
    {
      bavQt->getViewer()->m_time = false;

      pre_actionHG(bavQt->getViewer(), 'M');

      updateActiveChildTitle();
    }
  }
}

inline void
FEVV::SimpleWindow::on_actionAbout_MEPP_Help_triggered()
{
  QMessageBox::about(
      this,
      tr("About MEPP2 / Help"),
      tr("<b>MEPP2</b><br>"
         "<br>"
         "3D MEsh Processing Platform<br>"
         "Copyright (c) 2016-2019 University of Lyon and CNRS (France)<br>"
         "<br>"
         "LIRIS M2DISCO / MEPP-team<br>"
         "<br>"
         "<b>GitHub: <a href=\"https://github.com/MEPP-team/MEPP2\">see online "
         "repository</a></b><br>"
         "<b>Developer documentation: <a "
         "href=\"http://liris.cnrs.fr/mepp/doc/nightly/\">see online "
         "help</a></b><br>"
         "<br>"
         "-<br>"
         "<br>"
         "<b>Keys: </b><br>"
         "<br>"
         "Open in a new viewer -> <b>shift + OPEN</b><br>"
         "<br>"
         "Viewer -> <b>SELECT</b>         mesh : <b>shift</b><br>"
         "Viewer -> <b>TRANSLATE</b>      mesh : <b>T</b> (Translation "
         "draggers must be shown, see toolbar)<br>"
         "Viewer -> <b>ROTATE</b>         mesh : <b>R</b> (Rotation draggers "
         "must be shown, see toolbar)<br>"
         "<br>"
         "Viewer -> <b>OSG</b>            instrumentation : <b>S</b><br>"
         "<br>"
         "<b>*</b> on MacOS, the <b>ctrl</b> key is replaced by the "
         "<b>command</b> "
         "key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<br>"));
}

inline void
FEVV::SimpleWindow::onGrab()
{
#if(FEVV_USE_QT5)
  /*qApp->*/ QApplication::primaryScreen()->grabWindow(0).save(
      "screenshot.png");
#endif
}


inline
void
FEVV::SimpleWindow::actionHG(FEVV::SimpleViewer *viewer,
                             char t,
                             char t2)
{
  if(!viewer)
    return;

  if(t == 'A') // 'A'xis
  {
    viewer->gizmo->setNodeMask(viewer->m_ShowAxis ? 0xffffffff : 0x0);
  }
  else if(t == 'G') // 'G'rid
  {
    viewer->grid->setNodeMask(viewer->m_ShowGrid ? 0xffffffff : 0x0);
  }
  else if(t == '-') // TIME-
  {
    // time
    if(viewer->m_time)
    {
      std::vector< osg::Geode * > geodes = viewer->getGeodes();

      if(viewer->current_i_time >= 1)
      {
        if(t2 == '-')                 // TIME--
          viewer->current_i_time = 0; // FIRST ONE
        else
          viewer->current_i_time--; // PREVIOUS ONE

        std::vector< osg::Geode * > geodesSelected =
            viewer->getSelectedGeodes();
        for(unsigned i = 0; i < geodesSelected.size(); i++)
          viewer->setNodeSelected(geodesSelected[i], false);

        viewer->setNodeSelected(geodes[viewer->current_i_time], true);
        update(true);

        // std::cout << "current_i_time--" << std::endl;
      }
    }
    // time
  }
  else if(t == '+') // TIME+
  {
    // time
    if(viewer->m_time)
    {
      std::vector< osg::Geode * > geodes = viewer->getGeodes();

      if(geodes.size() >= 1 && viewer->current_i_time < (geodes.size() - 1))
      {
        if(t2 == '+')                              // TIME++
          viewer->current_i_time = viewer->i_time; // LAST ONE
        else
          viewer->current_i_time++; // NEXT ONE

        std::vector< osg::Geode * > geodesSelected =
            viewer->getSelectedGeodes();
        for(unsigned i = 0; i < geodesSelected.size(); i++)
          viewer->setNodeSelected(geodesSelected[i], false);

        viewer->setNodeSelected(geodes[viewer->current_i_time], true);
        update(true);

        // std::cout << "current_i_time++" << std::endl;
      }
    }
    // time
  }
  else if(t == 'M') // switch 'M'ODE - SPACE/TIME
  {
    // if (viewer->m_space_time) // BUT already protected in
    // 'on_actionChange_viewer_mode_triggered'
    {
      std::vector< osg::Geode * > geodes = viewer->getGeodes();

      int iSel = -1;
      bool force = false;

      for(unsigned i = 0; i < geodes.size(); i++)
      {
        if(viewer->isNodeSelected(geodes[i]))
        {
          iSel = i; // CURRENT ONE
          break;
        }
      }
      if(iSel == -1) // really IMPORTANT if MDIchild is active but no mesh
                     // selected in meshes list
      {
        if(viewer->current_i_time == -1)
          iSel = viewer->i_time; // LAST ONE
        else
          iSel = viewer->current_i_time; // SAVE ONE

        force = true;
      }

      std::vector< osg::Geode * > geodesSelected = viewer->getSelectedGeodes();
      for(unsigned i = 0; i < geodesSelected.size(); i++)
        viewer->setNodeSelected(geodesSelected[i], false);

      viewer->setNodeSelected(geodes[iSel], true);

      if((force) || (!viewer->isNodeSelected(geodes[iSel])))
        update(true);

      // ---

      if(viewer->m_time) // time
      {
        // std::cout << "entering in TIME MODE" << std::endl;

        viewer->current_i_time = iSel;
      }
      else // space
      {
        // std::cout << "entering in SPACE MODE" << std::endl;

        std::vector< osg::Group * > draggers1 = viewer->getDraggers1();
        std::vector< osg::Group * > draggers2 = viewer->getDraggers2();

        for(unsigned i = 0; i < geodes.size(); i++)
          geodes[i]->setNodeMask(0xffffffff);
        for(unsigned i = 0; i < draggers1.size(); i++)
          draggers1[i]->setNodeMask(viewer->m_ShowTranslateDragger ? 0xffffffff
                                                                   : 0x0);
        for(unsigned i = 0; i < draggers2.size(); i++)
          draggers2[i]->setNodeMask(viewer->m_ShowRotateDragger ? 0xffffffff
                                                                : 0x0);

        // actionHG(viewer, 'D', '_'); // 19/03/19 - test for calling a
        // re-'D'raw (because of 'hide/not hide' meshes in SPACE mode)
      }
    }
  }
  else if(t == 'T') // 'T'ranslate dragger
  {
    std::vector< osg::Group * > draggers1 = viewer->getDraggers1();

    if(viewer->m_time && viewer->current_i_time != -1)
    {
      for(unsigned i = 0; i < draggers1.size(); i++)
        draggers1[i]->setNodeMask(0x0);

      if(draggers1.size())
        draggers1[viewer->current_i_time]->setNodeMask(
            viewer->m_ShowTranslateDragger ? 0xffffffff : 0x0);
    }
    else
    {
      for(unsigned i = 0; i < draggers1.size(); i++)
        draggers1[i]->setNodeMask(viewer->m_ShowTranslateDragger ? 0xffffffff
                                                                 : 0x0);
    }
  }
  else if(t == 'R') // 'R'otate dragger
  {
    std::vector< osg::Group * > draggers2 = viewer->getDraggers2();

    if(viewer->m_time && viewer->current_i_time != -1)
    {
      for(unsigned i = 0; i < draggers2.size(); i++)
        draggers2[i]->setNodeMask(0x0);

      if(draggers2.size())
        draggers2[viewer->current_i_time]->setNodeMask(
            viewer->m_ShowRotateDragger ? 0xffffffff : 0x0);
    }
    else
    {
      for(unsigned i = 0; i < draggers2.size(); i++)
        draggers2[i]->setNodeMask(viewer->m_ShowRotateDragger ? 0xffffffff
                                                              : 0x0);
    }
  }
  else if(t == 'D') // re-'D'raw
  {
    FEVV::MixedMeshesVector meshes = viewer->getMeshes();
    std::vector< std::string > meshes_names = viewer->getMeshesNames();
    std::vector< FEVV::PMapsContainer * > properties_maps =
        viewer->get_properties_maps();

    for(unsigned i = 0; i < meshes.size(); i++)
    {
#ifdef FEVV_USE_CGAL
      if(meshes[i].second == "POLYHEDRON")
      {
        auto mesh_ptr = static_cast< FEVV::MeshPolyhedron* >(meshes[i].first);
        draw_or_redraw_mesh(mesh_ptr,
                            properties_maps[i],
                            viewer,
                            true,
                            false,
                            meshes_names[i],
                            false);
      }
      if(meshes[i].second == "SURFACEMESH")
      {
        auto mesh_ptr = static_cast< FEVV::MeshSurface* >(meshes[i].first);
        draw_or_redraw_mesh(mesh_ptr,
                            properties_maps[i],
                            viewer,
                            true,
                            false,
                            meshes_names[i],
                            false);
      }
      if(meshes[i].second == "LCC")
      {
        auto mesh_ptr = static_cast< FEVV::MeshLCC* >(meshes[i].first);
        draw_or_redraw_mesh(mesh_ptr,
                            properties_maps[i],
                            viewer,
                            true,
                            false,
                            meshes_names[i],
                            false);
      }
      if(meshes[i].second == "CGALPOINTSET")
      {
        auto mesh_ptr = static_cast< FEVV::CGALPointSet* >(meshes[i].first);
        draw_or_redraw_mesh(mesh_ptr,
                            properties_maps[i],
                            viewer,
                            true,
                            false,
                            meshes_names[i],
                            false);
      }
#endif //FEVV_USE_CGAL

#ifdef FEVV_USE_OPENMESH
      if(meshes[i].second == "OPENMESH")
      {
        auto mesh_ptr = static_cast< FEVV::MeshOpenMesh* >(meshes[i].first);
        draw_or_redraw_mesh(mesh_ptr,
                            properties_maps[i],
                            viewer,
                            true,
                            false,
                            meshes_names[i],
                            false);
      }
#endif //FEVV_USE_OPENMESH

#ifdef FEVV_USE_AIF
      if(meshes[i].second == "AIF")
      {
        auto mesh_ptr = static_cast< FEVV::MeshAIF* >(meshes[i].first);
        draw_or_redraw_mesh(mesh_ptr,
                            properties_maps[i],
                            viewer,
                            true,
                            false,
                            meshes_names[i],
                            false);
      }
#endif //FEVV_USE_AIF

#ifdef FEVV_USE_PCL
      if(meshes[i].second == "PCLPOINTCLOUD")
      {
        auto mesh_ptr = static_cast< FEVV::PCLPointCloud* >(meshes[i].first);
        draw_or_redraw_mesh(mesh_ptr,
                            properties_maps[i],
                            viewer,
                            true,
                            false,
                            meshes_names[i],
                            false);
      }
#endif //FEVV_USE_PCL
    }
  }
}


inline void
FEVV::SimpleWindow::pre_actionHG(Adapter::Viewer *base_viewer, char t, char t2)
{
  auto viewer = dynamic_cast< SimpleViewer * >(base_viewer);
  actionHG(viewer, t, t2);
}


inline void
FEVV::SimpleWindow::on_actionRender_Point_triggered()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

    bavQt->getViewer()->m_RenderMethod = RenderMethod::RENDER_POINTS;
    bavQt->getViewer()->m_space_time_changeColorMode = false;
    pre_actionHG(bavQt->getViewer());
    bavQt->getViewer()->m_space_time_changeColorMode = true;
  }
}

inline void
FEVV::SimpleWindow::on_actionRender_Line_triggered()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

    bavQt->getViewer()->m_RenderMethod = RenderMethod::RENDER_LINES;
    bavQt->getViewer()->m_space_time_changeColorMode = false;
    pre_actionHG(bavQt->getViewer());
    bavQt->getViewer()->m_space_time_changeColorMode = true;
  }
}

inline void
FEVV::SimpleWindow::on_actionRender_Fill_triggered()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

    bavQt->getViewer()->m_RenderMethod = RenderMethod::RENDER_FILL;
    bavQt->getViewer()->m_space_time_changeColorMode = false;
    pre_actionHG(bavQt->getViewer());
    bavQt->getViewer()->m_space_time_changeColorMode = true;
  }
}

inline void
FEVV::SimpleWindow::on_actionSuperimpose_Vertices_triggered()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

    bavQt->getViewer()->m_RenderSuperimposedVertices =
        !bavQt->getViewer()->m_RenderSuperimposedVertices;
    bavQt->getViewer()->m_space_time_changeColorMode = false;
    pre_actionHG(bavQt->getViewer());
    bavQt->getViewer()->m_space_time_changeColorMode = true;
  }
}

inline void
FEVV::SimpleWindow::on_actionSuperimpose_Vertices_bigger_triggered()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

    bavQt->getViewer()->m_RenderSuperimposedVertices_Big =
        !bavQt->getViewer()->m_RenderSuperimposedVertices_Big;
    bavQt->getViewer()->m_space_time_changeColorMode = false;
    pre_actionHG(bavQt->getViewer());
    bavQt->getViewer()->m_space_time_changeColorMode = true;
  }
}

inline void
FEVV::SimpleWindow::on_actionSuperimpose_Edges_triggered()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

    bavQt->getViewer()->m_RenderSuperimposedEdges =
        !bavQt->getViewer()->m_RenderSuperimposedEdges;
    bavQt->getViewer()->m_space_time_changeColorMode = false;
    pre_actionHG(bavQt->getViewer());
    bavQt->getViewer()->m_space_time_changeColorMode = true;
  }
}

inline void
FEVV::SimpleWindow::on_actionVertex_Color_triggered()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

    bavQt->getViewer()->m_UseVertexColor = true;
    bavQt->getViewer()->m_UseFaceColor = false;
    bavQt->getViewer()->m_UseTexture = false;
    pre_actionHG(bavQt->getViewer());
  }
}

inline void
FEVV::SimpleWindow::on_actionFace_Color_triggered()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

    bavQt->getViewer()->m_UseVertexColor = false;
    bavQt->getViewer()->m_UseFaceColor = true;
    bavQt->getViewer()->m_UseTexture = false;
    pre_actionHG(bavQt->getViewer());
  }
}

inline void
FEVV::SimpleWindow::on_actionTexture_Mode_triggered()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

    bavQt->getViewer()->m_UseVertexColor = false;
    bavQt->getViewer()->m_UseFaceColor = false;
    bavQt->getViewer()->m_UseTexture = true;
    pre_actionHG(bavQt->getViewer());
  }
}

inline void
FEVV::SimpleWindow::on_actionLighting_triggered()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

    bavQt->getViewer()->m_Lighting = !(bavQt->getViewer()->m_Lighting);
    bavQt->getViewer()->m_space_time_changeColorMode = false;
    pre_actionHG(bavQt->getViewer());
    bavQt->getViewer()->m_space_time_changeColorMode = true;
  }
}

inline void
FEVV::SimpleWindow::on_actionSmoothFlat_Shading_triggered()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

    bavQt->getViewer()->m_SmoothFlat_Shading =
        !bavQt->getViewer()->m_SmoothFlat_Shading;
    bavQt->getViewer()->m_space_time_changeColorMode = false;
    pre_actionHG(bavQt->getViewer());
    bavQt->getViewer()->m_space_time_changeColorMode = true;
  }
}

inline void
FEVV::SimpleWindow::on_actionRender_Mode_triggered()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

#if defined(__APPLE__)
    QMessageBox::information(
        this,
        "",
        QObject::tr("Only legacy rendering is available under OS X."));
#else
    bavQt->getViewer()->m_RenderMode = static_cast< RenderMode >(
        (static_cast< size_t >(bavQt->getViewer()->m_RenderMode) + 1) % 3);
    bavQt->getViewer()->m_space_time_changeColorMode = false;
    // bool SAVE_m_recreateOSGobj_if_redraw =
    // bavQt->getViewer()->m_recreateOSGobj_if_redraw;
    // bavQt->getViewer()->m_recreateOSGobj_if_redraw = true; // NEW
    pre_actionHG(bavQt->getViewer());
    // bavQt->getViewer()->m_recreateOSGobj_if_redraw =
    // SAVE_m_recreateOSGobj_if_redraw; // NEW
    bavQt->getViewer()->m_space_time_changeColorMode = true;

    updateActiveChildTitle();
#endif
  }
}


inline
void
FEVV::SimpleWindow::centerHG(FEVV::SimpleViewer *viewer)
{
  if(!viewer)
    return;

  FEVV::MixedMeshesVector meshes = viewer->getSelectedMeshes();

  for(unsigned i = 0; i < meshes.size(); i++)
  {
    viewer->centerMesh(meshes[i].first);
  }
}


inline void
FEVV::SimpleWindow::on_actionShow_Entire_Mesh_triggered()
{
  for(unsigned i = 0; i < adapters.size(); i++)
  {
    if(adapters[i]->isSelected())
    {
      BaseAdapterVisu *bav = adapters[i];
      auto viewer = dynamic_cast< SimpleViewer * >(bav->getViewer());
      centerHG(viewer);

    }
  }
}


inline void
FEVV::SimpleWindow::on_actionShow_Axis_triggered()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

    bavQt->getViewer()->m_ShowAxis = !(bavQt->getViewer()->m_ShowAxis);
    pre_actionHG(bavQt->getViewer(), 'A');
  }
}

inline void
FEVV::SimpleWindow::on_actionShow_Grid_triggered()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

    bavQt->getViewer()->m_ShowGrid = !(bavQt->getViewer()->m_ShowGrid);
    pre_actionHG(bavQt->getViewer(), 'G');
  }
}

inline void
FEVV::SimpleWindow::on_actionShow_Vertex_Normals_triggered()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

    bavQt->getViewer()->m_Show_Vertex_Normals =
        !bavQt->getViewer()->m_Show_Vertex_Normals;
    bavQt->getViewer()->m_space_time_changeColorMode = false;
    pre_actionHG(bavQt->getViewer());
    bavQt->getViewer()->m_space_time_changeColorMode = true;
  }
}

#if 0 //TODO-elo-rm-?-ask_MTO
template< typename HalfedgeGraph >
void
FEVV::SimpleWindow::showSelectedHG(FEVV::SimpleViewer *viewer)
{
  if(!viewer)
    return;
}
#endif


#if 0 //TODO-elo-rm-?-ask_MTO
inline void
FEVV::SimpleWindow::on_actionShow_Selected_triggered()
{
  for(unsigned i = 0; i < adapters.size(); i++)
  {
    if(adapters[i]->isSelected())
    {
      BaseAdapterVisu *bav = adapters[i];
      bav->getViewer()->m_ShowSelected = !(bav->getViewer()->m_ShowSelected);

#ifdef FEVV_USE_CGAL
      if(dynamic_cast< SimpleViewer< FEVV::MeshPolyhedron > * >(
             bav->getViewer()))
        showSelectedHG(dynamic_cast< SimpleViewer< FEVV::MeshPolyhedron > * >(
            bav->getViewer()));
      else if(dynamic_cast< SimpleViewer< FEVV::MeshSurface > * >(
                  bav->getViewer()))
        showSelectedHG(dynamic_cast< SimpleViewer< FEVV::MeshSurface > * >(
            bav->getViewer()));
      else if(dynamic_cast< SimpleViewer< FEVV::MeshLCC > * >(bav->getViewer()))
        showSelectedHG(
            dynamic_cast< SimpleViewer< FEVV::MeshLCC > * >(bav->getViewer()));
#endif
#ifdef FEVV_USE_OPENMESH
      if(dynamic_cast< SimpleViewer< FEVV::MeshOpenMesh > * >(bav->getViewer()))
        showSelectedHG(dynamic_cast< SimpleViewer< FEVV::MeshOpenMesh > * >(
            bav->getViewer()));
#endif
#ifdef FEVV_USE_AIF
      if(dynamic_cast< SimpleViewer< FEVV::MeshAIF > * >(bav->getViewer()))
        showSelectedHG(
            dynamic_cast< SimpleViewer< FEVV::MeshAIF > * >(bav->getViewer()));
#endif
    }
  }
}
#endif


inline void
FEVV::SimpleWindow::on_actionShow_Translation_Draggers_triggered()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

    bavQt->getViewer()->m_ShowTranslateDragger =
        !(bavQt->getViewer()->m_ShowTranslateDragger);
    pre_actionHG(bavQt->getViewer(), 'T');
  }
}

inline void
FEVV::SimpleWindow::on_actionShow_Rotation_Draggers_triggered()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

    bavQt->getViewer()->m_ShowRotateDragger =
        !(bavQt->getViewer()->m_ShowRotateDragger);
    pre_actionHG(bavQt->getViewer(), 'R');
  }
}

inline void
FEVV::SimpleWindow::on_actionDyn_First_triggered()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

    if(bavQt->getViewer()->m_space_time)
      if(bavQt->getViewer()->m_time)
        pre_actionHG(bavQt->getViewer(), '-', '-');
  }
}

inline void
FEVV::SimpleWindow::on_actionDyn_Previous_triggered()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

    if(bavQt->getViewer()->m_space_time)
      if(bavQt->getViewer()->m_time)
        pre_actionHG(bavQt->getViewer(), '-');
  }
}

inline void
FEVV::SimpleWindow::on_actionDyn_Next_triggered()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

    if(bavQt->getViewer()->m_space_time)
      if(bavQt->getViewer()->m_time)
        pre_actionHG(bavQt->getViewer(), '+');
  }
}

inline void
FEVV::SimpleWindow::on_actionDyn_Last_triggered()
{
  if(activeMdiChild())
  {
    BaseAdapterVisuQt *bavQt =
        dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());

    if(bavQt->getViewer()->m_space_time)
      if(bavQt->getViewer()->m_time)
        pre_actionHG(bavQt->getViewer(), '+', '+');
  }
}

inline void
FEVV::SimpleWindow::onAddBall()
{
  if(!Assert::check(isValid(),
                    "is not valid (see init() or attach()). Leaving...",
                    "SimpleWindow::onAddBall"))
  {
    return;
  }

  if(adapters.size() > 0)
  {
    int x = (std::rand() % 20) - 10;
    int y = (std::rand() % 20) - 10;
    int z = (std::rand() % 20) - 10;

    std::cout << "Adding a ball in {" << x << "," << y << "," << z << "}"
              << std::endl;

    BaseAdapterVisuQt *bavQt = dynamic_cast< BaseAdapterVisuQt * >(adapters[0]);
    static_cast< BaseViewerOSG * >(adapters[0]->getViewer())
        ->addGroup(Debug::createBall(
            x,
            y,
            z,
            1,
            Color::Amethyst(),
            "Ball [" + bavQt->windowTitle().toStdString() + std::string("]")));
  }
}

inline void
FEVV::SimpleWindow::onItemChanged(QListWidgetItem *_item,
                                  QListWidgetItem *_item_old)
{
  using Viewer = BaseViewerOSG;

  /// removing old selected item
  if(_item_old != 0)
  {
    QVariant variant = _item_old->data(Qt::UserRole);
    Viewer::DataModel data = variant.value< Viewer::DataModel >();
    if(data.type == Helpers::DataType::MODEL && data.node != nullptr)
    {
      data.viewer->setNodeSelected(data.node, false);
      data.viewer->setSelected(false);
      for_each(adapters.begin(), adapters.end(), [&data](Adapter *w) {
        if(w->getViewer() == data.viewer)
        {
          w->setSelected(false);
        }
      });
    }
  }

  /// adding new selected item
  if(_item != 0)
  {
    QVariant variant = _item->data(Qt::UserRole);
    Viewer::DataModel data = variant.value< Viewer::DataModel >();
    if(data.type == Helpers::DataType::MODEL && data.node != nullptr)
    {
      data.viewer->setNodeSelected(data.node, true);
      data.viewer->setSelected(true);
      for_each(adapters.begin(), adapters.end(), [&data](Adapter *w) {
        if(w->getViewer() == data.viewer)
        {
          w->setSelected(true);

          // QMdiSubWindow FOCUS
          SimpleWindow *sw =
              dynamic_cast< SimpleWindow * >(w->getViewer()->getWindow());
          QMdiArea *_mdiArea = sw->mdiArea;
          if(_mdiArea)
          {
            BaseAdapterVisuQt *bavQt = dynamic_cast< BaseAdapterVisuQt * >(w);
            QList< QMdiSubWindow * > listMdiSubW = _mdiArea->subWindowList();

            for(int MdiSubW = 0; MdiSubW < listMdiSubW.size(); ++MdiSubW)
            {
              BaseAdapterVisuQt *bavQtMdiSubW =
                  dynamic_cast< BaseAdapterVisuQt * >(
                      listMdiSubW.at(MdiSubW)->widget());
              if(bavQtMdiSubW->getViewer() == bavQt->getViewer())
              {
                _mdiArea->setActiveSubWindow(listMdiSubW.at(MdiSubW));
                // std::cout << "onItemChanged(): " <<
                // listMdiSubW.at(MdiSubW)->windowTitle().toStdString() <<
                // std::endl;
                break;
              }
            }
          }
          // QMdiSubWindow FOCUS
        }
      });
    }
  }
}

inline std::vector< FEVV::SimpleWindow::Adapter * >
FEVV::SimpleWindow::getSelectedAdapters()
{
  std::vector< Adapter * > result;
  for_each(adapters.begin(), adapters.end(), [&result](Adapter *w) {
    if(w->isSelected())
    {
      result.push_back(w);
    }
  });

  return result;
}

inline std::vector< FEVV::SimpleWindow::Adapter::Viewer * >
FEVV::SimpleWindow::getSelectedViewers()
{
  std::vector< Adapter::Viewer * > result;
  for_each(adapters.begin(), adapters.end(), [&result](Adapter *w) {
    if(w->isSelected())
    {
      result.push_back(w->getViewer());
    }
  });

  return result;
}

inline
std::string
FEVV::SimpleWindow::chooseDatastructureMsgBox(void)
{
  QMessageBox msgbox;
  msgbox.setWindowTitle("Datastructure");
  msgbox.setText("Choose a datastructure to store the mesh(es):");
  msgbox.setIcon(QMessageBox::Question);

  // here the Role is used to order the buttons
#ifdef FEVV_USE_CGAL
  QPushButton *polyhedron_button =
      msgbox.addButton("Polyhedron_3", QMessageBox::ResetRole);
  QPushButton *surfacemesh_button =
      msgbox.addButton("Surface_mesh", QMessageBox::ResetRole);
  QPushButton *lcc_button =
      msgbox.addButton("LCC", QMessageBox::ResetRole);
  QPushButton *cgalpointset_button =
      msgbox.addButton("CGALPointSet", QMessageBox::ResetRole);
#endif

#ifdef FEVV_USE_OPENMESH
  QPushButton *openmesh_button =
      msgbox.addButton("OpenMesh", QMessageBox::ResetRole);
#endif

#ifdef FEVV_USE_AIF
  QPushButton *aif_button = msgbox.addButton("AIF", QMessageBox::ResetRole);
#endif

#ifdef FEVV_USE_PCL
  QPushButton *pcl_button =
      msgbox.addButton("PCLPointCloud", QMessageBox::ResetRole);
#endif

  QPushButton *abortButton = msgbox.addButton(QMessageBox::Cancel);

  msgbox.exec();
  std::string choice("NONE");

#ifdef FEVV_USE_CGAL
  if(msgbox.clickedButton() == polyhedron_button)
  {
    choice = "POLYHEDRON";
  }
  else if(msgbox.clickedButton() == surfacemesh_button)
  {
    choice = "SURFACEMESH";
  }
  else if(msgbox.clickedButton() == lcc_button)
  {
    choice = "LCC";
  }
  else if(msgbox.clickedButton() == cgalpointset_button)
  {
    choice = "CGALPOINTSET";
  }
#endif

#ifdef FEVV_USE_OPENMESH
  if(msgbox.clickedButton() == openmesh_button)
  {
    choice = "OPENMESH";
  }
#endif

#ifdef FEVV_USE_AIF
  if(msgbox.clickedButton() == aif_button)
  {
    choice = "AIF";
  }
#endif

#ifdef FEVV_USE_PCL
  if(msgbox.clickedButton() == pcl_button)
  {
    choice = "PCLPOINTCLOUD";
  }
#endif

  return choice;
}


inline
FEVV::SimpleViewer *
FEVV::SimpleWindow::createNewViewer(void)
{
  FEVV::SimpleAdapterVisu *adapter;
  adapter= new FEVV::SimpleAdapterVisu();
  adapter->setWindowTitle(QObject::tr("<aid: %1>")
                          .arg((qlonglong)adapter, 0, 16));
  adapter->setWindowIcon(QIcon(":/logo/resources/MEPP.png"));
  FEVV::SimpleViewer *viewer = new FEVV::SimpleViewer();
  adapter->attach(viewer);
  adapter->init();
  viewer->init();
  viewer->setSelected(true);
  viewer->m_space_time = true;
  enableSpaceTimeMenus();
  this->attach(adapter, USE_MDI);
  adapter->show();
  if(this->getMdiArea())
    this->getMdiArea()->tileSubWindows();

  return viewer;
}
