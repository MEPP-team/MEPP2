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
#endif

inline FEVV::SimpleWindow::SimpleWindow(QWidget *_parent,
                                        Qt::WindowFlags _flags)
    : BaseWindowQt(_parent, _flags)
{
  ui.setupUi(this);
}

inline FEVV::SimpleWindow::~SimpleWindow()
{
  if(mdiArea != nullptr)
  {
    delete mdiArea;
  }
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
  if(!isValid())
  {
    QMessageBox::information(this, "",
		QObject::tr("Please, first <b>open a mesh</b><br>"
					"<br>"
					"or if necessary with some plugins <b>use an empty child window</b><br>"
					"(see 'Open -> EMPTY mode' <b>key</b> in 'About MEPP2 / Help' menu)."));

    return;
  }

  for(auto funct : stackPlugins)
  {
    for_each(adapters.begin(), adapters.end(), [&funct](Adapter *w) {
      if(w->isSelected())
      {
        w->apply(funct.second);
      }
    });
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
FEVV::SimpleWindow::load_mesh(const std::string &mesh_filename,
                              MeshT *&mesh,
                              FEVV::PMapsContainer &pmaps_bag)
{
  statusBar()->showMessage(QObject::tr("Open mesh file...") /*, 2000*/);

  mesh = new MeshT;
  FEVV::Filters::read_mesh(mesh_filename, *mesh, pmaps_bag, open_only_pts_mode);

  statusBar()->showMessage(QObject::tr("") /*, 2000*/);
}

template< typename MeshT >
void
FEVV::SimpleWindow::empty_mesh(MeshT *&mesh, FEVV::PMapsContainer &pmaps_bag)
{
  mesh = new MeshT;
}

template< typename MeshT >
void
FEVV::SimpleWindow::draw_or_redraw_mesh(
    /*const */ MeshT *mesh,
    /*const */ FEVV::PMapsContainer *pmaps_bag,
    FEVV::SimpleViewer< MeshT > *viewer,
    bool _redraw,
    bool _recomputeNT_if_redraw,
    std::string _mesh_filename,
    float _step)
{
  // QApplication::setOverrideCursor(Qt::BusyCursor);
  viewer->draw_or_redraw_mesh(mesh, pmaps_bag, _redraw, _recomputeNT_if_redraw, _mesh_filename, _step);
  // QApplication::restoreOverrideCursor();
}

inline void
FEVV::SimpleWindow::on_actionNew_triggered()
{
  QMessageBox::information(this, "", QObject::tr("New - Not yet implemented."));
}

inline void
FEVV::SimpleWindow::on_actionOpen_Polyhedron_3_SPACE_TIME()
{
#ifdef FEVV_USE_CGAL
  QString validExtensions =
      "OBJ/OFF Files (*.obj *.off);;OBJ Files (*.obj);;OFF files (*.off);;COFF "
      "files (*.coff);;PLY files (*.ply);;MSH files (*.msh)";
  QString validVtkExtensions =
      "VTK Files (*.vtk);;VTP files (*.vtp);;VTU files (*.vtu)";

  QString allExtensions = validExtensions
#ifdef FEVV_USE_VTK
                          + ";;" + validVtkExtensions
#endif
#ifdef FEVV_USE_FBX
                          + ";;FBX files (*.fbx)"
#endif
                          + "";

  QString suffix;
  QFileDialog::Option options = (QFileDialog::Option)0;

#if defined(__linux__) || defined(__APPLE__)
  options = QFileDialog::DontUseNativeDialog; // PB under LINUX !?
#endif

  QStringList files =
      QFileDialog::getOpenFileNames(this,
                                    "Open Polyhedron_3 (SPACE/TIME)",
                                    /*openLocation*/ QDir::currentPath(),
                                    allExtensions,
                                    &suffix,
                                    options);

  ///// Polyhedron
  FEVV::SimpleAdapterVisu< FEVV::MeshPolyhedron >
      *adapterPolyhedron; // if (!USE_MDI) delete
  FEVV::SimpleViewer< FEVV::MeshPolyhedron >
      *viewerPolyhedron; // already destroy by the adapter destructor

  if(files.size() > 0)
  {
    adapterPolyhedron = new FEVV::SimpleAdapterVisu< FEVV::MeshPolyhedron >();
    // adapterPolyhedron->setWindowTitle(QObject::tr("Polyhedron_3"));
    adapterPolyhedron->setWindowTitle(
        QObject::tr("<Polyhedron_3 - aid: %1>")
            .arg((qlonglong)adapterPolyhedron, 0, 16));
    adapterPolyhedron->setWindowIcon(QIcon(":/logo/resources/MEPP.png"));
    viewerPolyhedron = new FEVV::SimpleViewer< FEVV::MeshPolyhedron >();
    adapterPolyhedron->attach(viewerPolyhedron);
    adapterPolyhedron->init();
    viewerPolyhedron->init();
    viewerPolyhedron->m_space_time = true;
    enableSpaceTimeMenus();
    this->attach(adapterPolyhedron, USE_MDI);
  }

  QStringList::Iterator it = files.begin();
  int m = 0;
  while(it != files.end())
  {
    FEVV::MeshPolyhedron *mPolyhedron =
        nullptr; // already destroy by the viewer destructor
    FEVV::PMapsContainer *p_polyhedron_pmaps_bag =
        new FEVV::PMapsContainer; // already destroy by the viewer destructor

    QApplication::setOverrideCursor(Qt::WaitCursor);
    load_mesh((*it).toStdString(), mPolyhedron, *p_polyhedron_pmaps_bag);
    QApplication::restoreOverrideCursor();

    draw_or_redraw_mesh(
        mPolyhedron,
        p_polyhedron_pmaps_bag,
        viewerPolyhedron,
        false,
        false,
        FEVV::FileUtils::get_file_full_name((*it).toStdString()),
        m * STEP_SPACE);

    ++it;
    ++m;
  }

  if(files.size() > 0)
  {
    adapterPolyhedron->show();
    updateActiveChildTitle();
    if(this->getMdiArea())
      this->getMdiArea()->tileSubWindows();
  }
#else
  QMessageBox::information(
      this,
      "",
      QObject::tr("actionOpen_Polyhedron_3 : FEVV_USE_CGAL is OFF !"));
#endif
}

#ifdef FEVV_USE_CGAL
inline void
FEVV::SimpleWindow::on_actionOpen_Polyhedron_3_ADD_SPACE_TIME(
    FEVV::SimpleViewer< FEVV::MeshPolyhedron > *viewerPolyhedron)
{
#ifdef FEVV_USE_CGAL
  QString validExtensions =
      "OBJ/OFF Files (*.obj *.off);;OBJ Files (*.obj);;OFF files (*.off);;COFF "
      "files (*.coff);;PLY files (*.ply);;MSH files (*.msh)";
  QString validVtkExtensions =
      "VTK Files (*.vtk);;VTP files (*.vtp);;VTU files (*.vtu)";

  QString allExtensions = validExtensions
#ifdef FEVV_USE_VTK
                          + ";;" + validVtkExtensions
#endif
#ifdef FEVV_USE_FBX
                          + ";;FBX files (*.fbx)"
#endif
                          + "";

  QString suffix;
  QFileDialog::Option options = (QFileDialog::Option)0;

#if defined(__linux__) || defined(__APPLE__)
  options = QFileDialog::DontUseNativeDialog; // PB under LINUX !?
#endif

  QStringList files =
      QFileDialog::getOpenFileNames(this,
                                    "Open Polyhedron_3 (ADD SPACE/TIME)",
                                    /*openLocation*/ QDir::currentPath(),
                                    allExtensions,
                                    &suffix,
                                    options);

  if(files.size() > 0)
  {
    viewerPolyhedron->m_space_time = true;
    enableSpaceTimeMenus();
    updateActiveChildTitle();
  }

  QStringList::Iterator it = files.begin();
  int m = 0;
  while(it != files.end())
  {
    FEVV::MeshPolyhedron *mPolyhedron =
        nullptr; // already destroy by the viewer destructor
    FEVV::PMapsContainer *p_polyhedron_pmaps_bag =
        new FEVV::PMapsContainer; // already destroy by the viewer destructor

    QApplication::setOverrideCursor(Qt::WaitCursor);
    load_mesh((*it).toStdString(), mPolyhedron, *p_polyhedron_pmaps_bag);
    QApplication::restoreOverrideCursor();

    draw_or_redraw_mesh(
        mPolyhedron,
        p_polyhedron_pmaps_bag,
        viewerPolyhedron,
        false,
        false,
        FEVV::FileUtils::get_file_full_name((*it).toStdString()),
        m * STEP_SPACE);

    ++it;
    ++m;
  }
#else
  QMessageBox::information(
      this,
      "",
      QObject::tr("actionOpen_Polyhedron_3 : FEVV_USE_CGAL is OFF !"));
#endif
}
#endif

inline void
FEVV::SimpleWindow::on_actionOpen_Polyhedron_3_NORMAL(bool empty)
{
#ifdef FEVV_USE_CGAL
  QString validExtensions =
      "OBJ/OFF Files (*.obj *.off);;OBJ Files (*.obj);;OFF files (*.off);;COFF "
      "files (*.coff);;PLY files (*.ply);;MSH files (*.msh)";
  QString validVtkExtensions =
      "VTK Files (*.vtk);;VTP files (*.vtp);;VTU files (*.vtu)";

  QString allExtensions = validExtensions
#ifdef FEVV_USE_VTK
                          + ";;" + validVtkExtensions
#endif
#ifdef FEVV_USE_FBX
                          + ";;FBX files (*.fbx)"
#endif
                          + "";

  QString suffix;
  QFileDialog::Option options = (QFileDialog::Option)0;

#if defined(__linux__) || defined(__APPLE__)
  options = QFileDialog::DontUseNativeDialog; // PB under LINUX !?
#endif

  QStringList files;

  if(!empty)
  {
    files = QFileDialog::getOpenFileNames(this,
                                          "Open Polyhedron_3",
                                          /*openLocation*/ QDir::currentPath(),
                                          allExtensions,
                                          &suffix,
                                          options);
  }

  QStringList::Iterator it = files.begin();
  while(empty || it != files.end())
  {
    ///// Polyhedron
    FEVV::SimpleAdapterVisu< FEVV::MeshPolyhedron >
        *adapterPolyhedron; // if (!USE_MDI) delete
    FEVV::SimpleViewer< FEVV::MeshPolyhedron >
        *viewerPolyhedron; // already destroy by the adapter destructor

    FEVV::MeshPolyhedron *mPolyhedron =
        nullptr; // already destroy by the viewer destructor
    FEVV::PMapsContainer *p_polyhedron_pmaps_bag =
        new FEVV::PMapsContainer; // already destroy by the viewer destructor

    adapterPolyhedron = new FEVV::SimpleAdapterVisu< FEVV::MeshPolyhedron >();
    // adapterPolyhedron->setWindowTitle(QObject::tr("Polyhedron_3"));
    adapterPolyhedron->setWindowTitle(
        QObject::tr("<Polyhedron_3 - aid: %1>")
            .arg((qlonglong)adapterPolyhedron, 0, 16));
    adapterPolyhedron->setWindowIcon(QIcon(":/logo/resources/MEPP.png"));
    viewerPolyhedron = new FEVV::SimpleViewer< FEVV::MeshPolyhedron >();
    adapterPolyhedron->attach(viewerPolyhedron);
    adapterPolyhedron->init();
    viewerPolyhedron->init();
    this->attach(adapterPolyhedron, USE_MDI);

    std::string mesh_filename = "";

    QApplication::setOverrideCursor(Qt::WaitCursor);
    if(empty)
    {
      empty_mesh(mPolyhedron, *p_polyhedron_pmaps_bag);
      mesh_filename = "empty";
    }
    else
    {
      load_mesh((*it).toStdString(), mPolyhedron, *p_polyhedron_pmaps_bag);
      mesh_filename = FEVV::FileUtils::get_file_full_name((*it).toStdString());
    }
    QApplication::restoreOverrideCursor();

    draw_or_redraw_mesh(mPolyhedron,
                        p_polyhedron_pmaps_bag,
                        viewerPolyhedron,
                        false,
                        false,
                        mesh_filename);

    adapterPolyhedron->show();
    updateActiveChildTitle();
    if(this->getMdiArea())
      this->getMdiArea()->tileSubWindows();

    if(empty)
      empty = false;
    else
      ++it;
  }
#else
  QMessageBox::information(
      this,
      "",
      QObject::tr("actionOpen_Polyhedron_3 : FEVV_USE_CGAL is OFF !"));
#endif
}

inline void
FEVV::SimpleWindow::on_actionOpen_Polyhedron_3_triggered()
{
  if(QApplication::keyboardModifiers().testFlag(Qt::AltModifier)) // for 'only_pts' mode
    open_only_pts_mode = true;
  else
    open_only_pts_mode = false;

#ifdef __APPLE__
  // on macOS, the ControlModifier value corresponds to the Command keys on the
  // keyboard, and the MetaModifier value corresponds to the Control keys
  if(QApplication::keyboardModifiers().testFlag(Qt::ShiftModifier) &&
     QApplication::keyboardModifiers().testFlag(
         /*Qt::MetaModifier*/ Qt::ControlModifier))
#else
  if(QApplication::keyboardModifiers().testFlag(Qt::ShiftModifier) &&
     QApplication::keyboardModifiers().testFlag(Qt::ControlModifier))
#endif
  {
    if(activeMdiChild())
    {
      BaseAdapterVisuQt *bavQt =
          dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());
      Adapter::Viewer *viewer = bavQt->getViewer();

#ifdef FEVV_USE_CGAL
      if(dynamic_cast< SimpleViewer< FEVV::MeshPolyhedron > * >(viewer))
        on_actionOpen_Polyhedron_3_ADD_SPACE_TIME(
            dynamic_cast< SimpleViewer< FEVV::MeshPolyhedron > * >(viewer));
#else
      QMessageBox::information(
          this,
          "",
          QObject::tr("actionOpen_Polyhedron_3 : FEVV_USE_CGAL is OFF !"));
#endif
    }
  }
  else if(QApplication::keyboardModifiers().testFlag(
              Qt::ShiftModifier)) // QApplication::queryKeyboardModifiers()
    on_actionOpen_Polyhedron_3_SPACE_TIME();
#ifdef __APPLE__
  // on macOS, the ControlModifier value corresponds to the Command keys on the
  // keyboard, and the MetaModifier value corresponds to the Control keys
  else if(QApplication::keyboardModifiers().testFlag(
              /*Qt::MetaModifier*/ Qt::ControlModifier))
#else
  else if(QApplication::keyboardModifiers().testFlag(Qt::ControlModifier))
#endif
    on_actionOpen_Polyhedron_3_NORMAL(true); // empty mesh
  else
    on_actionOpen_Polyhedron_3_NORMAL();
}

inline void
FEVV::SimpleWindow::on_actionOpen_Surface_mesh_SPACE_TIME()
{
#ifdef FEVV_USE_CGAL
  QString validExtensions =
      "OBJ/OFF Files (*.obj *.off);;OBJ Files (*.obj);;OFF files (*.off);;COFF "
      "files (*.coff);;PLY files (*.ply);;MSH files (*.msh)";
  QString validVtkExtensions =
      "VTK Files (*.vtk);;VTP files (*.vtp);;VTU files (*.vtu)";

  QString allExtensions = validExtensions
#ifdef FEVV_USE_VTK
                          + ";;" + validVtkExtensions
#endif
#ifdef FEVV_USE_FBX
                          + ";;FBX files (*.fbx)"
#endif
                          + "";

  QString suffix;
  QFileDialog::Option options = (QFileDialog::Option)0;

#if defined(__linux__) || defined(__APPLE__)
  options = QFileDialog::DontUseNativeDialog; // PB under LINUX !?
#endif

  QStringList files =
      QFileDialog::getOpenFileNames(this,
                                    "Open Surface_mesh (SPACE/TIME)",
                                    /*openLocation*/ QDir::currentPath(),
                                    allExtensions,
                                    &suffix,
                                    options);

  ///// Surface_mesh
  FEVV::SimpleAdapterVisu< FEVV::MeshSurface >
      *adapterSurface; // if (!USE_MDI) delete
  FEVV::SimpleViewer< FEVV::MeshSurface >
      *viewerSurface; // already destroy by the adapter destructor

  if(files.size() > 0)
  {
    adapterSurface = new FEVV::SimpleAdapterVisu< FEVV::MeshSurface >();
    // adapterSurface->setWindowTitle(QObject::tr("Surface_mesh"));
    adapterSurface->setWindowTitle(QObject::tr("<Surface_mesh - aid: %1>")
                                       .arg((qlonglong)adapterSurface, 0, 16));
    adapterSurface->setWindowIcon(QIcon(":/logo/resources/MEPP.png"));
    viewerSurface = new FEVV::SimpleViewer< FEVV::MeshSurface >();
    adapterSurface->attach(viewerSurface);
    adapterSurface->init();
    viewerSurface->init();
    viewerSurface->m_space_time = true;
    enableSpaceTimeMenus();
    this->attach(adapterSurface, USE_MDI);
  }

  QStringList::Iterator it = files.begin();
  int m = 0;
  while(it != files.end())
  {
    FEVV::MeshSurface *mSurface =
        nullptr; // already destroy by the viewer destructor
    FEVV::PMapsContainer *p_surfacemesh_pmaps_bag =
        new FEVV::PMapsContainer; // already destroy by the viewer destructor

    QApplication::setOverrideCursor(Qt::WaitCursor);
    load_mesh((*it).toStdString(), mSurface, *p_surfacemesh_pmaps_bag);
    QApplication::restoreOverrideCursor();

    draw_or_redraw_mesh(
        mSurface,
        p_surfacemesh_pmaps_bag,
        viewerSurface,
        false,
        false,
        FEVV::FileUtils::get_file_full_name((*it).toStdString()),
        m * STEP_SPACE);

    ++it;
    ++m;
  }

  if(files.size() > 0)
  {
    adapterSurface->show();
    updateActiveChildTitle();
    if(this->getMdiArea())
      this->getMdiArea()->tileSubWindows();
  }
#else
  QMessageBox::information(
      this,
      "",
      QObject::tr("actionOpen_Surface_mesh : FEVV_USE_CGAL is OFF !"));
#endif
}

#ifdef FEVV_USE_CGAL
inline void
FEVV::SimpleWindow::on_actionOpen_Surface_mesh_ADD_SPACE_TIME(
    FEVV::SimpleViewer< FEVV::MeshSurface > *viewerSurface)
{
#ifdef FEVV_USE_CGAL
  QString validExtensions =
      "OBJ/OFF Files (*.obj *.off);;OBJ Files (*.obj);;OFF files (*.off);;COFF "
      "files (*.coff);;PLY files (*.ply);;MSH files (*.msh)";
  QString validVtkExtensions =
      "VTK Files (*.vtk);;VTP files (*.vtp);;VTU files (*.vtu)";

  QString allExtensions = validExtensions
#ifdef FEVV_USE_VTK
                          + ";;" + validVtkExtensions
#endif
#ifdef FEVV_USE_FBX
                          + ";;FBX files (*.fbx)"
#endif
                          + "";

  QString suffix;
  QFileDialog::Option options = (QFileDialog::Option)0;

#if defined(__linux__) || defined(__APPLE__)
  options = QFileDialog::DontUseNativeDialog; // PB under LINUX !?
#endif

  QStringList files =
      QFileDialog::getOpenFileNames(this,
                                    "Open Surface_mesh (ADD SPACE/TIME)",
                                    /*openLocation*/ QDir::currentPath(),
                                    allExtensions,
                                    &suffix,
                                    options);

  if(files.size() > 0)
  {
    viewerSurface->m_space_time = true;
    enableSpaceTimeMenus();
    updateActiveChildTitle();
  }

  QStringList::Iterator it = files.begin();
  int m = 0;
  while(it != files.end())
  {
    FEVV::MeshSurface *mSurface =
        nullptr; // already destroy by the viewer destructor
    FEVV::PMapsContainer *p_surfacemesh_pmaps_bag =
        new FEVV::PMapsContainer; // already destroy by the viewer destructor

    QApplication::setOverrideCursor(Qt::WaitCursor);
    load_mesh((*it).toStdString(), mSurface, *p_surfacemesh_pmaps_bag);
    QApplication::restoreOverrideCursor();

    draw_or_redraw_mesh(
        mSurface,
        p_surfacemesh_pmaps_bag,
        viewerSurface,
        false,
        false,
        FEVV::FileUtils::get_file_full_name((*it).toStdString()),
        m * STEP_SPACE);

    ++it;
    ++m;
  }
#else
  QMessageBox::information(
      this,
      "",
      QObject::tr("actionOpen_Surface_mesh : FEVV_USE_CGAL is OFF !"));
#endif
}
#endif

inline void
FEVV::SimpleWindow::on_actionOpen_Surface_mesh_NORMAL(bool empty)
{
#ifdef FEVV_USE_CGAL
  QString validExtensions =
      "OBJ/OFF Files (*.obj *.off);;OBJ Files (*.obj);;OFF files (*.off);;COFF "
      "files (*.coff);;PLY files (*.ply);;MSH files (*.msh)";
  QString validVtkExtensions =
      "VTK Files (*.vtk);;VTP files (*.vtp);;VTU files (*.vtu)";

  QString allExtensions = validExtensions
#ifdef FEVV_USE_VTK
                          + ";;" + validVtkExtensions
#endif
#ifdef FEVV_USE_FBX
                          + ";;FBX files (*.fbx)"
#endif
                          + "";

  QString suffix;
  QFileDialog::Option options = (QFileDialog::Option)0;

#if defined(__linux__) || defined(__APPLE__)
  options = QFileDialog::DontUseNativeDialog; // PB under LINUX !?
#endif

  QStringList files;

  if(!empty)
  {
    files = QFileDialog::getOpenFileNames(this,
                                          "Open Surface_mesh",
                                          /*openLocation*/ QDir::currentPath(),
                                          allExtensions,
                                          &suffix,
                                          options);
  }

  QStringList::Iterator it = files.begin();
  while(empty || it != files.end())
  {
    ///// Surface_mesh
    FEVV::SimpleAdapterVisu< FEVV::MeshSurface >
        *adapterSurface; // if (!USE_MDI) delete
    FEVV::SimpleViewer< FEVV::MeshSurface >
        *viewerSurface; // already destroy by the adapter destructor

    FEVV::MeshSurface *mSurface =
        nullptr; // already destroy by the viewer destructor
    FEVV::PMapsContainer *p_surfacemesh_pmaps_bag =
        new FEVV::PMapsContainer; // already destroy by the viewer destructor

    adapterSurface = new FEVV::SimpleAdapterVisu< FEVV::MeshSurface >();
    // adapterSurface->setWindowTitle(QObject::tr("Surface_mesh"));
    adapterSurface->setWindowTitle(QObject::tr("<Surface_mesh - aid: %1>")
                                       .arg((qlonglong)adapterSurface, 0, 16));
    adapterSurface->setWindowIcon(QIcon(":/logo/resources/MEPP.png"));
    viewerSurface = new FEVV::SimpleViewer< FEVV::MeshSurface >();
    adapterSurface->attach(viewerSurface);
    adapterSurface->init();
    viewerSurface->init();
    this->attach(adapterSurface, USE_MDI);

    std::string mesh_filename = "";

    QApplication::setOverrideCursor(Qt::WaitCursor);
    if(empty)
    {
      empty_mesh(mSurface, *p_surfacemesh_pmaps_bag);
      mesh_filename = "empty";
    }
    else
    {
      load_mesh((*it).toStdString(), mSurface, *p_surfacemesh_pmaps_bag);
      mesh_filename = FEVV::FileUtils::get_file_full_name((*it).toStdString());
    }
    QApplication::restoreOverrideCursor();

    draw_or_redraw_mesh(
        mSurface, p_surfacemesh_pmaps_bag, viewerSurface, false, false, mesh_filename);

    adapterSurface->show();
    updateActiveChildTitle();
    if(this->getMdiArea())
      this->getMdiArea()->tileSubWindows();

    if(empty)
      empty = false;
    else
      ++it;
  }
#else
  QMessageBox::information(
      this,
      "",
      QObject::tr("actionOpen_Surface_mesh : FEVV_USE_CGAL is OFF !"));
#endif
}

inline void
FEVV::SimpleWindow::on_actionOpen_Surface_mesh_triggered()
{
  if(QApplication::keyboardModifiers().testFlag(Qt::AltModifier)) // for 'only_pts' mode
    open_only_pts_mode = true;
  else
    open_only_pts_mode = false;

#ifdef __APPLE__
  // on macOS, the ControlModifier value corresponds to the Command keys on the
  // keyboard, and the MetaModifier value corresponds to the Control keys
  if(QApplication::keyboardModifiers().testFlag(Qt::ShiftModifier) &&
     QApplication::keyboardModifiers().testFlag(
         /*Qt::MetaModifier*/ Qt::ControlModifier))
#else
  if(QApplication::keyboardModifiers().testFlag(Qt::ShiftModifier) &&
     QApplication::keyboardModifiers().testFlag(Qt::ControlModifier))
#endif
  {
    if(activeMdiChild())
    {
      BaseAdapterVisuQt *bavQt =
          dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());
      Adapter::Viewer *viewer = bavQt->getViewer();

#ifdef FEVV_USE_CGAL
      if(dynamic_cast< SimpleViewer< FEVV::MeshSurface > * >(viewer))
        on_actionOpen_Surface_mesh_ADD_SPACE_TIME(
            dynamic_cast< SimpleViewer< FEVV::MeshSurface > * >(viewer));
#else
      QMessageBox::information(
          this,
          "",
          QObject::tr("actionOpen_Surface_mesh : FEVV_USE_CGAL is OFF !"));
#endif
    }
  }
  else if(QApplication::keyboardModifiers().testFlag(
              Qt::ShiftModifier)) // QApplication::queryKeyboardModifiers()
    on_actionOpen_Surface_mesh_SPACE_TIME();
#ifdef __APPLE__
  // on macOS, the ControlModifier value corresponds to the Command keys on the
  // keyboard, and the MetaModifier value corresponds to the Control keys
  else if(QApplication::keyboardModifiers().testFlag(
              /*Qt::MetaModifier*/ Qt::ControlModifier))
#else
  else if(QApplication::keyboardModifiers().testFlag(Qt::ControlModifier))
#endif
    on_actionOpen_Surface_mesh_NORMAL(true); // empty mesh
  else
    on_actionOpen_Surface_mesh_NORMAL();
}

inline void
FEVV::SimpleWindow::on_actionOpen_Linear_cell_complex_SPACE_TIME()
{
#ifdef FEVV_USE_CGAL

  QString validExtensions =
      "OBJ/OFF Files (*.obj *.off);;OBJ Files (*.obj);;OFF files (*.off);;COFF "
      "files (*.coff);;PLY files (*.ply);;MSH files (*.msh)";
  QString validVtkExtensions =
      "VTK Files (*.vtk);;VTP files (*.vtp);;VTU files (*.vtu)";

  QString allExtensions = validExtensions
#ifdef FEVV_USE_VTK
                          + ";;" + validVtkExtensions
#endif
#ifdef FEVV_USE_FBX
                          + ";;FBX files (*.fbx)"
#endif
                          + "";

  QString suffix;
  QFileDialog::Option options = (QFileDialog::Option)0;

#if defined(__linux__) || defined(__APPLE__)
  options = QFileDialog::DontUseNativeDialog; // PB under LINUX !?
#endif

  QStringList files =
      QFileDialog::getOpenFileNames(this,
                                    "Open Linear_cell_complex (SPACE/TIME)",
                                    /*openLocation*/ QDir::currentPath(),
                                    allExtensions,
                                    &suffix,
                                    options);

  ///// Linear_cell_complex
  FEVV::SimpleAdapterVisu< FEVV::MeshLCC > *adapterLCC; // if (!USE_MDI) delete
  FEVV::SimpleViewer< FEVV::MeshLCC >
      *viewerLCC; // already destroy by the adapter destructor

  if(files.size() > 0)
  {
    adapterLCC = new FEVV::SimpleAdapterVisu< FEVV::MeshLCC >();
    // adapterLCC->setWindowTitle(QObject::tr("Linear_cell_complex"));
    adapterLCC->setWindowTitle(QObject::tr("<Linear_cell_complex - aid: %1>")
                                   .arg((qlonglong)adapterLCC, 0, 16));
    adapterLCC->setWindowIcon(QIcon(":/logo/resources/MEPP.png"));
    viewerLCC = new FEVV::SimpleViewer< FEVV::MeshLCC >();
    adapterLCC->attach(viewerLCC);
    adapterLCC->init();
    viewerLCC->init();
    viewerLCC->m_space_time = true;
    enableSpaceTimeMenus();
    this->attach(adapterLCC, USE_MDI);
  }

  QStringList::Iterator it = files.begin();
  int m = 0;
  while(it != files.end())
  {
    FEVV::MeshLCC *mLCC = nullptr; // already destroy by the viewer destructor
    FEVV::PMapsContainer *p_lcc_pmaps_bag =
        new FEVV::PMapsContainer; // already destroy by the viewer destructor

    QApplication::setOverrideCursor(Qt::WaitCursor);
    load_mesh((*it).toStdString(), mLCC, *p_lcc_pmaps_bag);
    QApplication::restoreOverrideCursor();

    draw_or_redraw_mesh(
        mLCC,
        p_lcc_pmaps_bag,
        viewerLCC,
        false,
        false,
        FEVV::FileUtils::get_file_full_name((*it).toStdString()),
        m * STEP_SPACE);

    ++it;
    ++m;
  }

  if(files.size() > 0)
  {
    adapterLCC->show();
    updateActiveChildTitle();
    if(this->getMdiArea())
      this->getMdiArea()->tileSubWindows();
  }
#else
  QMessageBox::information(
      this,
      "",
      QObject::tr("actionOpen_Linear_cell_complex : FEVV_USE_CGAL is OFF !"));
#endif
}

#ifdef FEVV_USE_CGAL
inline void
FEVV::SimpleWindow::on_actionOpen_Linear_cell_complex_ADD_SPACE_TIME(
    FEVV::SimpleViewer< FEVV::MeshLCC > *viewerLCC)
{
#ifdef FEVV_USE_CGAL

  QString validExtensions =
      "OBJ/OFF Files (*.obj *.off);;OBJ Files (*.obj);;OFF files (*.off);;COFF "
      "files (*.coff);;PLY files (*.ply);;MSH files (*.msh)";
  QString validVtkExtensions =
      "VTK Files (*.vtk);;VTP files (*.vtp);;VTU files (*.vtu)";

  QString allExtensions = validExtensions
#ifdef FEVV_USE_VTK
                          + ";;" + validVtkExtensions
#endif
#ifdef FEVV_USE_FBX
                          + ";;FBX files (*.fbx)"
#endif
                          + "";

  QString suffix;
  QFileDialog::Option options = (QFileDialog::Option)0;

#if defined(__linux__) || defined(__APPLE__)
  options = QFileDialog::DontUseNativeDialog; // PB under LINUX !?
#endif

  QStringList files =
      QFileDialog::getOpenFileNames(this,
                                    "Open Linear_cell_complex (ADD SPACE/TIME)",
                                    /*openLocation*/ QDir::currentPath(),
                                    allExtensions,
                                    &suffix,
                                    options);

  if(files.size() > 0)
  {
    viewerLCC->m_space_time = true;
    enableSpaceTimeMenus();
    updateActiveChildTitle();
  }

  QStringList::Iterator it = files.begin();
  int m = 0;
  while(it != files.end())
  {
    FEVV::MeshLCC *mLCC = nullptr; // already destroy by the viewer destructor
    FEVV::PMapsContainer *p_lcc_pmaps_bag =
        new FEVV::PMapsContainer; // already destroy by the viewer destructor

    QApplication::setOverrideCursor(Qt::WaitCursor);
    load_mesh((*it).toStdString(), mLCC, *p_lcc_pmaps_bag);
    QApplication::restoreOverrideCursor();

    draw_or_redraw_mesh(
        mLCC,
        p_lcc_pmaps_bag,
        viewerLCC,
        false,
        false,
        FEVV::FileUtils::get_file_full_name((*it).toStdString()),
        m * STEP_SPACE);

    ++it;
    ++m;
  }
#else
  QMessageBox::information(
      this,
      "",
      QObject::tr("actionOpen_Linear_cell_complex : FEVV_USE_CGAL is OFF !"));
#endif
}
#endif

inline void
FEVV::SimpleWindow::on_actionOpen_Linear_cell_complex_NORMAL(bool empty)
{
#ifdef FEVV_USE_CGAL

  QString validExtensions =
      "OBJ/OFF Files (*.obj *.off);;OBJ Files (*.obj);;OFF files (*.off);;COFF "
      "files (*.coff);;PLY files (*.ply);;MSH files (*.msh)";
  QString validVtkExtensions =
      "VTK Files (*.vtk);;VTP files (*.vtp);;VTU files (*.vtu)";

  QString allExtensions = validExtensions
#ifdef FEVV_USE_VTK
                          + ";;" + validVtkExtensions
#endif
#ifdef FEVV_USE_FBX
                          + ";;FBX files (*.fbx)"
#endif
                          + "";

  QString suffix;
  QFileDialog::Option options = (QFileDialog::Option)0;

#if defined(__linux__) || defined(__APPLE__)
  options = QFileDialog::DontUseNativeDialog; // PB under LINUX !?
#endif

  QStringList files;

  if(!empty)
  {
    files = QFileDialog::getOpenFileNames(this,
                                          "Open Linear_cell_complex",
                                          /*openLocation*/ QDir::currentPath(),
                                          allExtensions,
                                          &suffix,
                                          options);
  }

  QStringList::Iterator it = files.begin();
  while(empty || it != files.end())
  {
    ///// Linear_cell_complex
    FEVV::SimpleAdapterVisu< FEVV::MeshLCC >
        *adapterLCC; // if (!USE_MDI) delete
    FEVV::SimpleViewer< FEVV::MeshLCC >
        *viewerLCC; // already destroy by the adapter destructor

    FEVV::MeshLCC *mLCC = nullptr; // already destroy by the viewer destructor
    FEVV::PMapsContainer *p_lcc_pmaps_bag =
        new FEVV::PMapsContainer; // already destroy by the viewer destructor

    adapterLCC = new FEVV::SimpleAdapterVisu< FEVV::MeshLCC >();
    // adapterLCC->setWindowTitle(QObject::tr("Linear_cell_complex"));
    adapterLCC->setWindowTitle(QObject::tr("<Linear_cell_complex - aid: %1>")
                                   .arg((qlonglong)adapterLCC, 0, 16));
    adapterLCC->setWindowIcon(QIcon(":/logo/resources/MEPP.png"));
    viewerLCC = new FEVV::SimpleViewer< FEVV::MeshLCC >();
    adapterLCC->attach(viewerLCC);
    adapterLCC->init();
    viewerLCC->init();
    this->attach(adapterLCC, USE_MDI);

    std::string mesh_filename = "";

    QApplication::setOverrideCursor(Qt::WaitCursor);
    if(empty)
    {
      empty_mesh(mLCC, *p_lcc_pmaps_bag);
      mesh_filename = "empty";
    }
    else
    {
      load_mesh((*it).toStdString(), mLCC, *p_lcc_pmaps_bag);
      mesh_filename = FEVV::FileUtils::get_file_full_name((*it).toStdString());
    }
    QApplication::restoreOverrideCursor();

    draw_or_redraw_mesh(mLCC, p_lcc_pmaps_bag, viewerLCC, false, false, mesh_filename);

    adapterLCC->show();
    updateActiveChildTitle();
    if(this->getMdiArea())
      this->getMdiArea()->tileSubWindows();

    if(empty)
      empty = false;
    else
      ++it;
  }
#else
  QMessageBox::information(
      this,
      "",
      QObject::tr("actionOpen_Linear_cell_complex : FEVV_USE_CGAL is OFF !"));
#endif
}

inline void
FEVV::SimpleWindow::on_actionOpen_Linear_cell_complex_triggered()
{
  if(QApplication::keyboardModifiers().testFlag(Qt::AltModifier)) // for 'only_pts' mode
    open_only_pts_mode = true;
  else
    open_only_pts_mode = false;

#ifdef __APPLE__
  // on macOS, the ControlModifier value corresponds to the Command keys on the
  // keyboard, and the MetaModifier value corresponds to the Control keys
  if(QApplication::keyboardModifiers().testFlag(Qt::ShiftModifier) &&
     QApplication::keyboardModifiers().testFlag(
         /*Qt::MetaModifier*/ Qt::ControlModifier))
#else
  if(QApplication::keyboardModifiers().testFlag(Qt::ShiftModifier) &&
     QApplication::keyboardModifiers().testFlag(Qt::ControlModifier))
#endif
  {
    if(activeMdiChild())
    {
      BaseAdapterVisuQt *bavQt =
          dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());
      Adapter::Viewer *viewer = bavQt->getViewer();

#ifdef FEVV_USE_CGAL
      if(dynamic_cast< SimpleViewer< FEVV::MeshLCC > * >(viewer))
        on_actionOpen_Linear_cell_complex_ADD_SPACE_TIME(
            dynamic_cast< SimpleViewer< FEVV::MeshLCC > * >(viewer));
#else
      QMessageBox::information(
          this,
          "",
          QObject::tr(
              "actionOpen_Linear_cell_complex : FEVV_USE_CGAL is OFF !"));
#endif
    }
  }
  else if(QApplication::keyboardModifiers().testFlag(
              Qt::ShiftModifier)) // QApplication::queryKeyboardModifiers()
    on_actionOpen_Linear_cell_complex_SPACE_TIME();
#ifdef __APPLE__
  // on macOS, the ControlModifier value corresponds to the Command keys on the
  // keyboard, and the MetaModifier value corresponds to the Control keys
  else if(QApplication::keyboardModifiers().testFlag(
              /*Qt::MetaModifier*/ Qt::ControlModifier))
#else
  else if(QApplication::keyboardModifiers().testFlag(Qt::ControlModifier))
#endif
    on_actionOpen_Linear_cell_complex_NORMAL(true); // empty mesh
  else
    on_actionOpen_Linear_cell_complex_NORMAL();
}

inline void
FEVV::SimpleWindow::on_actionOpen_OpenMesh_SPACE_TIME()
{
#ifdef FEVV_USE_OPENMESH
  QString validExtensions =
      "OBJ/OFF Files (*.obj *.off);;OBJ Files (*.obj);;OFF files (*.off);;COFF "
      "files (*.coff);;PLY files (*.ply);;MSH files (*.msh)";
  QString validVtkExtensions =
      "VTK Files (*.vtk);;VTP files (*.vtp);;VTU files (*.vtu)";

  QString allExtensions = validExtensions
#ifdef FEVV_USE_VTK
                          + ";;" + validVtkExtensions
#endif
#ifdef FEVV_USE_FBX
                          + ";;FBX files (*.fbx)"
#endif
                          + "";

  QString suffix;
  QFileDialog::Option options = (QFileDialog::Option)0;

#if defined(__linux__) || defined(__APPLE__)
  options = QFileDialog::DontUseNativeDialog; // PB under LINUX !?
#endif

  QStringList files =
      QFileDialog::getOpenFileNames(this,
                                    "Open OpenMesh (SPACE/TIME)",
                                    /*openLocation*/ QDir::currentPath(),
                                    allExtensions,
                                    &suffix,
                                    options);

  ///// OpenMesh
  FEVV::SimpleAdapterVisu< FEVV::MeshOpenMesh >
      *adapterOpenMesh; // if (!USE_MDI) delete
  FEVV::SimpleViewer< FEVV::MeshOpenMesh >
      *viewerOpenMesh; // already destroy by the adapter destructor

  if(files.size() > 0)
  {
    adapterOpenMesh = new FEVV::SimpleAdapterVisu< FEVV::MeshOpenMesh >();
    // adapterOpenMesh->setWindowTitle(QObject::tr("OpenMesh"));
    adapterOpenMesh->setWindowTitle(
        QObject::tr("<OpenMesh - aid: %1>")
            .arg((qlonglong)adapterOpenMesh, 0, 16));
    adapterOpenMesh->setWindowIcon(QIcon(":/logo/resources/MEPP.png"));
    viewerOpenMesh = new FEVV::SimpleViewer< FEVV::MeshOpenMesh >();
    adapterOpenMesh->attach(viewerOpenMesh);
    adapterOpenMesh->init();
    viewerOpenMesh->init();
    viewerOpenMesh->m_space_time = true;
    enableSpaceTimeMenus();
    this->attach(adapterOpenMesh, USE_MDI);
  }

  QStringList::Iterator it = files.begin();
  int m = 0;
  while(it != files.end())
  {
    FEVV::MeshOpenMesh *mOpenMesh =
        nullptr; // already destroy by the viewer destructor
    FEVV::PMapsContainer *p_openmesh_pmaps_bag =
        new FEVV::PMapsContainer; // already destroy by the viewer destructor

    QApplication::setOverrideCursor(Qt::WaitCursor);
    load_mesh((*it).toStdString(), mOpenMesh, *p_openmesh_pmaps_bag);
    QApplication::restoreOverrideCursor();

    draw_or_redraw_mesh(
        mOpenMesh,
        p_openmesh_pmaps_bag,
        viewerOpenMesh,
        false,
        false,
        FEVV::FileUtils::get_file_full_name((*it).toStdString()),
        m * STEP_SPACE);

    ++it;
    ++m;
  }

  if(files.size() > 0)
  {
    adapterOpenMesh->show();
    updateActiveChildTitle();
    if(this->getMdiArea())
      this->getMdiArea()->tileSubWindows();
  }
#else
  QMessageBox::information(
      this, "", QObject::tr("actionOpen_OpenMesh : FEVV_USE_CGAL is OFF !"));
#endif
}

#ifdef FEVV_USE_OPENMESH
inline void
FEVV::SimpleWindow::on_actionOpen_OpenMesh_ADD_SPACE_TIME(
    FEVV::SimpleViewer< FEVV::MeshOpenMesh > *viewerOpenMesh)
{
#ifdef FEVV_USE_OPENMESH
  QString validExtensions =
      "OBJ/OFF Files (*.obj *.off);;OBJ Files (*.obj);;OFF files (*.off);;COFF "
      "files (*.coff);;PLY files (*.ply);;MSH files (*.msh)";
  QString validVtkExtensions =
      "VTK Files (*.vtk);;VTP files (*.vtp);;VTU files (*.vtu)";

  QString allExtensions = validExtensions
#ifdef FEVV_USE_VTK
                          + ";;" + validVtkExtensions
#endif
#ifdef FEVV_USE_FBX
                          + ";;FBX files (*.fbx)"
#endif
                          + "";

  QString suffix;
  QFileDialog::Option options = (QFileDialog::Option)0;

#if defined(__linux__) || defined(__APPLE__)
  options = QFileDialog::DontUseNativeDialog; // PB under LINUX !?
#endif

  QStringList files =
      QFileDialog::getOpenFileNames(this,
                                    "Open OpenMesh (ADD SPACE/TIME)",
                                    /*openLocation*/ QDir::currentPath(),
                                    allExtensions,
                                    &suffix,
                                    options);

  if(files.size() > 0)
  {
    viewerOpenMesh->m_space_time = true;
    enableSpaceTimeMenus();
    updateActiveChildTitle();
  }

  QStringList::Iterator it = files.begin();
  int m = 0;
  while(it != files.end())
  {
    FEVV::MeshOpenMesh *mOpenMesh =
        nullptr; // already destroy by the viewer destructor
    FEVV::PMapsContainer *p_openmesh_pmaps_bag =
        new FEVV::PMapsContainer; // already destroy by the viewer destructor

    QApplication::setOverrideCursor(Qt::WaitCursor);
    load_mesh((*it).toStdString(), mOpenMesh, *p_openmesh_pmaps_bag);
    QApplication::restoreOverrideCursor();

    draw_or_redraw_mesh(
        mOpenMesh,
        p_openmesh_pmaps_bag,
        viewerOpenMesh,
        false,
        false,
        FEVV::FileUtils::get_file_full_name((*it).toStdString()),
        m * STEP_SPACE);

    ++it;
    ++m;
  }
#else
  QMessageBox::information(
      this, "", QObject::tr("actionOpen_OpenMesh : FEVV_USE_CGAL is OFF !"));
#endif
}
#endif

inline void
FEVV::SimpleWindow::on_actionOpen_OpenMesh_NORMAL(bool empty)
{
#ifdef FEVV_USE_OPENMESH
  QString validExtensions =
      "OBJ/OFF Files (*.obj *.off);;OBJ Files (*.obj);;OFF files (*.off);;COFF "
      "files (*.coff);;PLY files (*.ply);;MSH files (*.msh)";
  QString validVtkExtensions =
      "VTK Files (*.vtk);;VTP files (*.vtp);;VTU files (*.vtu)";

  QString allExtensions = validExtensions
#ifdef FEVV_USE_VTK
                          + ";;" + validVtkExtensions
#endif
#ifdef FEVV_USE_FBX
                          + ";;FBX files (*.fbx)"
#endif
                          + "";

  QString suffix;
  QFileDialog::Option options = (QFileDialog::Option)0;

#if defined(__linux__) || defined(__APPLE__)
  options = QFileDialog::DontUseNativeDialog; // PB under LINUX !?
#endif

  QStringList files;

  if(!empty)
  {
    files = QFileDialog::getOpenFileNames(this,
                                          "Open OpenMesh",
                                          /*openLocation*/ QDir::currentPath(),
                                          allExtensions,
                                          &suffix,
                                          options);
  }

  QStringList::Iterator it = files.begin();
  while(empty || it != files.end())
  {
    ///// OpenMesh
    FEVV::SimpleAdapterVisu< FEVV::MeshOpenMesh >
        *adapterOpenMesh; // if (!USE_MDI) delete
    FEVV::SimpleViewer< FEVV::MeshOpenMesh >
        *viewerOpenMesh; // already destroy by the adapter destructor

    FEVV::MeshOpenMesh *mOpenMesh =
        nullptr; // already destroy by the viewer destructor
    FEVV::PMapsContainer *p_openmesh_pmaps_bag =
        new FEVV::PMapsContainer; // already destroy by the viewer destructor

    adapterOpenMesh = new FEVV::SimpleAdapterVisu< FEVV::MeshOpenMesh >();
    // adapterOpenMesh->setWindowTitle(QObject::tr("OpenMesh"));
    adapterOpenMesh->setWindowTitle(
        QObject::tr("<OpenMesh - aid: %1>")
            .arg((qlonglong)adapterOpenMesh, 0, 16));
    adapterOpenMesh->setWindowIcon(QIcon(":/logo/resources/MEPP.png"));
    viewerOpenMesh = new FEVV::SimpleViewer< FEVV::MeshOpenMesh >();
    adapterOpenMesh->attach(viewerOpenMesh);
    adapterOpenMesh->init();
    viewerOpenMesh->init();
    this->attach(adapterOpenMesh, USE_MDI);

    std::string mesh_filename = "";

    QApplication::setOverrideCursor(Qt::WaitCursor);
    if(empty)
    {
      empty_mesh(mOpenMesh, *p_openmesh_pmaps_bag);
      mesh_filename = "empty";
    }
    else
    {
      load_mesh((*it).toStdString(), mOpenMesh, *p_openmesh_pmaps_bag);
      mesh_filename = FEVV::FileUtils::get_file_full_name((*it).toStdString());
    }
    QApplication::restoreOverrideCursor();

    draw_or_redraw_mesh(
        mOpenMesh, p_openmesh_pmaps_bag, viewerOpenMesh, false, false, mesh_filename);

    adapterOpenMesh->show();
    updateActiveChildTitle();
    if(this->getMdiArea())
      this->getMdiArea()->tileSubWindows();

    if(empty)
      empty = false;
    else
      ++it;
  }
#else
  QMessageBox::information(
      this,
      "",
      QObject::tr("actionOpen_OpenMesh : FEVV_USE_OPENMESH is OFF !"));
#endif
}

inline void
FEVV::SimpleWindow::on_actionOpen_OpenMesh_triggered()
{
  if(QApplication::keyboardModifiers().testFlag(Qt::AltModifier)) // for 'only_pts' mode
    open_only_pts_mode = true;
  else
    open_only_pts_mode = false;

#ifdef __APPLE__
  // on macOS, the ControlModifier value corresponds to the Command keys on the
  // keyboard, and the MetaModifier value corresponds to the Control keys
  if(QApplication::keyboardModifiers().testFlag(Qt::ShiftModifier) &&
     QApplication::keyboardModifiers().testFlag(
         /*Qt::MetaModifier*/ Qt::ControlModifier))
#else
  if(QApplication::keyboardModifiers().testFlag(Qt::ShiftModifier) &&
     QApplication::keyboardModifiers().testFlag(Qt::ControlModifier))
#endif
  {
    if(activeMdiChild())
    {
      BaseAdapterVisuQt *bavQt =
          dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());
      Adapter::Viewer *viewer = bavQt->getViewer();

#ifdef FEVV_USE_OPENMESH
      if(dynamic_cast< SimpleViewer< FEVV::MeshOpenMesh > * >(viewer))
        on_actionOpen_OpenMesh_ADD_SPACE_TIME(
            dynamic_cast< SimpleViewer< FEVV::MeshOpenMesh > * >(viewer));
#else
      QMessageBox::information(
          this,
          "",
          QObject::tr("actionOpen_OpenMesh : FEVV_USE_OPENMESH is OFF !"));
#endif
    }
  }
  else if(QApplication::keyboardModifiers().testFlag(
              Qt::ShiftModifier)) // QApplication::queryKeyboardModifiers()
    on_actionOpen_OpenMesh_SPACE_TIME();
#ifdef __APPLE__
  // on macOS, the ControlModifier value corresponds to the Command keys on the
  // keyboard, and the MetaModifier value corresponds to the Control keys
  else if(QApplication::keyboardModifiers().testFlag(
              /*Qt::MetaModifier*/ Qt::ControlModifier))
#else
  else if(QApplication::keyboardModifiers().testFlag(Qt::ControlModifier))
#endif
    on_actionOpen_OpenMesh_NORMAL(true); // empty mesh
  else
    on_actionOpen_OpenMesh_NORMAL();
}

inline void
FEVV::SimpleWindow::on_actionOpen_AIFMesh_SPACE_TIME()
{
#ifdef FEVV_USE_AIF
  QString validExtensions =
      "OBJ/OFF Files (*.obj *.off);;OBJ Files (*.obj);;OFF files (*.off);;COFF "
      "files (*.coff);;PLY files (*.ply);;MSH files (*.msh)";
  QString validVtkExtensions =
      "VTK Files (*.vtk);;VTP files (*.vtp);;VTU files (*.vtu)";

  QString allExtensions = validExtensions
#ifdef FEVV_USE_VTK
                          + ";;" + validVtkExtensions
#endif
#ifdef FEVV_USE_FBX
                          + ";;FBX files (*.fbx)"
#endif
                          + "";

  QString suffix;
  QFileDialog::Option options = (QFileDialog::Option)0;

#if defined(__linux__) || defined(__APPLE__)
  options = QFileDialog::DontUseNativeDialog; // PB under LINUX !?
#endif

  QStringList files =
      QFileDialog::getOpenFileNames(this,
                                    "Open AIFMesh (SPACE/TIME)",
                                    /*openLocation*/ QDir::currentPath(),
                                    allExtensions,
                                    &suffix,
                                    options);

  ///// AIFMesh
  FEVV::SimpleAdapterVisu< FEVV::MeshAIF > *adapterAIF; // if (!USE_MDI) delete
  FEVV::SimpleViewer< FEVV::MeshAIF >
      *viewerAIF; // already destroy by the adapter destructor

  if(files.size() > 0)
  {
    adapterAIF = new FEVV::SimpleAdapterVisu< FEVV::MeshAIF >();
    // adapterAIF->setWindowTitle(QObject::tr("AIFMesh"));
    adapterAIF->setWindowTitle(
        QObject::tr("<AIFMesh - aid: %1>").arg((qlonglong)adapterAIF, 0, 16));
    adapterAIF->setWindowIcon(QIcon(":/logo/resources/MEPP.png"));
    viewerAIF = new FEVV::SimpleViewer< FEVV::MeshAIF >();
    adapterAIF->attach(viewerAIF);
    adapterAIF->init();
    viewerAIF->init();
    viewerAIF->m_space_time = true;
    enableSpaceTimeMenus();
    this->attach(adapterAIF, USE_MDI);
  }

  QStringList::Iterator it = files.begin();
  int m = 0;
  while(it != files.end())
  {
    FEVV::MeshAIF *mAIF = nullptr; // already destroy by the viewer destructor
    FEVV::PMapsContainer *p_aif_pmaps_bag =
        new FEVV::PMapsContainer; // already destroy by the viewer destructor

    QApplication::setOverrideCursor(Qt::WaitCursor);
    load_mesh((*it).toStdString(), mAIF, *p_aif_pmaps_bag);
    QApplication::restoreOverrideCursor();

    draw_or_redraw_mesh(
        mAIF,
        p_aif_pmaps_bag,
        viewerAIF,
        false,
        false,
        FEVV::FileUtils::get_file_full_name((*it).toStdString()),
        m * STEP_SPACE);

    ++it;
    ++m;
  }

  if(files.size() > 0)
  {
    adapterAIF->show();
    updateActiveChildTitle();
    if(this->getMdiArea())
      this->getMdiArea()->tileSubWindows();
  }
#else
  QMessageBox::information(
      this, "", QObject::tr("actionOpen_AIFMesh : FEVV_USE_CGAL is OFF !"));
#endif
}

#ifdef FEVV_USE_AIF
inline void
FEVV::SimpleWindow::on_actionOpen_AIFMesh_ADD_SPACE_TIME(
    FEVV::SimpleViewer< FEVV::MeshAIF > *viewerAIF)
{
#ifdef FEVV_USE_AIF
  QString validExtensions =
      "OBJ/OFF Files (*.obj *.off);;OBJ Files (*.obj);;OFF files (*.off);;COFF "
      "files (*.coff);;PLY files (*.ply);;MSH files (*.msh)";
  QString validVtkExtensions =
      "VTK Files (*.vtk);;VTP files (*.vtp);;VTU files (*.vtu)";

  QString allExtensions = validExtensions
#ifdef FEVV_USE_VTK
                          + ";;" + validVtkExtensions
#endif
#ifdef FEVV_USE_FBX
                          + ";;FBX files (*.fbx)"
#endif
                          + "";

  QString suffix;
  QFileDialog::Option options = (QFileDialog::Option)0;

#if defined(__linux__) || defined(__APPLE__)
  options = QFileDialog::DontUseNativeDialog; // PB under LINUX !?
#endif

  QStringList files =
      QFileDialog::getOpenFileNames(this,
                                    "Open AIFMesh (ADD SPACE/TIME)",
                                    /*openLocation*/ QDir::currentPath(),
                                    allExtensions,
                                    &suffix,
                                    options);

  if(files.size() > 0)
  {
    viewerAIF->m_space_time = true;
    enableSpaceTimeMenus();
    updateActiveChildTitle();
  }

  QStringList::Iterator it = files.begin();
  int m = 0;
  while(it != files.end())
  {
    FEVV::MeshAIF *mAIF = nullptr; // already destroy by the viewer destructor
    FEVV::PMapsContainer *p_aif_pmaps_bag =
        new FEVV::PMapsContainer; // already destroy by the viewer destructor

    QApplication::setOverrideCursor(Qt::WaitCursor);
    load_mesh((*it).toStdString(), mAIF, *p_aif_pmaps_bag);
    QApplication::restoreOverrideCursor();

    draw_or_redraw_mesh(
        mAIF,
        p_aif_pmaps_bag,
        viewerAIF,
        false,
        false,
        FEVV::FileUtils::get_file_full_name((*it).toStdString()),
        m * STEP_SPACE);

    ++it;
    ++m;
  }
#else
  QMessageBox::information(
      this, "", QObject::tr("actionOpen_AIFMesh : FEVV_USE_CGAL is OFF !"));
#endif
}
#endif

inline void
FEVV::SimpleWindow::on_actionOpen_AIFMesh_NORMAL(bool empty)
{
#ifdef FEVV_USE_AIF
  QString validExtensions =
      "OBJ/OFF Files (*.obj *.off);;OBJ Files (*.obj);;OFF files (*.off);;COFF "
      "files (*.coff);;PLY files (*.ply);;MSH files (*.msh)";
  QString validVtkExtensions =
      "VTK Files (*.vtk);;VTP files (*.vtp);;VTU files (*.vtu)";

  QString allExtensions = validExtensions
#ifdef FEVV_USE_VTK
                          + ";;" + validVtkExtensions
#endif
#ifdef FEVV_USE_FBX
                          + ";;FBX files (*.fbx)"
#endif
                          + "";

  QString suffix;
  QFileDialog::Option options = (QFileDialog::Option)0;

#if defined(__linux__) || defined(__APPLE__)
  options = QFileDialog::DontUseNativeDialog; // PB under LINUX !?
#endif

  QStringList files;

  if(!empty)
  {
    files = QFileDialog::getOpenFileNames(this,
                                          "Open AIFMesh",
                                          /*openLocation*/ QDir::currentPath(),
                                          allExtensions,
                                          &suffix,
                                          options);
  }

  QStringList::Iterator it = files.begin();
  while(empty || it != files.end())
  {
    ///// AIFMesh
    FEVV::SimpleAdapterVisu< FEVV::MeshAIF >
        *adapterAIF; // if (!USE_MDI) delete
    FEVV::SimpleViewer< FEVV::MeshAIF >
        *viewerAIF; // already destroy by the adapter destructor

    FEVV::MeshAIF *mAIF = nullptr; // already destroy by the viewer destructor
    FEVV::PMapsContainer *p_aif_pmaps_bag =
        new FEVV::PMapsContainer; // already destroy by the viewer destructor

    adapterAIF = new FEVV::SimpleAdapterVisu< FEVV::MeshAIF >();
    // adapterAIF->setWindowTitle(QObject::tr("AIFMesh"));
    adapterAIF->setWindowTitle(
        QObject::tr("<AIFMesh - aid: %1>").arg((qlonglong)adapterAIF, 0, 16));
    adapterAIF->setWindowIcon(QIcon(":/logo/resources/MEPP.png"));
    viewerAIF = new FEVV::SimpleViewer< FEVV::MeshAIF >();
    adapterAIF->attach(viewerAIF);
    adapterAIF->init();
    viewerAIF->init();
    this->attach(adapterAIF, USE_MDI);

    std::string mesh_filename = "";

    QApplication::setOverrideCursor(Qt::WaitCursor);
    if(empty)
    {
      empty_mesh(mAIF, *p_aif_pmaps_bag);
      mesh_filename = "empty";
    }
    else
    {
      load_mesh((*it).toStdString(), mAIF, *p_aif_pmaps_bag);
      mesh_filename = FEVV::FileUtils::get_file_full_name((*it).toStdString());
    }
    QApplication::restoreOverrideCursor();

    draw_or_redraw_mesh(mAIF, p_aif_pmaps_bag, viewerAIF, false, false, mesh_filename);

    adapterAIF->show();
    updateActiveChildTitle();
    if(this->getMdiArea())
      this->getMdiArea()->tileSubWindows();

    if(empty)
      empty = false;
    else
      ++it;
  }
#else
  QMessageBox::information(
      this, "", QObject::tr("actionOpen_AIFMesh : FEVV_USE_AIF is OFF !"));
#endif
}

inline void
FEVV::SimpleWindow::on_actionOpen_AIFMesh_triggered()
{
  if(QApplication::keyboardModifiers().testFlag(Qt::AltModifier)) // for 'only_pts' mode
    open_only_pts_mode = true;
  else
    open_only_pts_mode = false;

#ifdef __APPLE__
  // on macOS, the ControlModifier value corresponds to the Command keys on the
  // keyboard, and the MetaModifier value corresponds to the Control keys
  if(QApplication::keyboardModifiers().testFlag(Qt::ShiftModifier) &&
     QApplication::keyboardModifiers().testFlag(
         /*Qt::MetaModifier*/ Qt::ControlModifier))
#else
  if(QApplication::keyboardModifiers().testFlag(Qt::ShiftModifier) &&
     QApplication::keyboardModifiers().testFlag(Qt::ControlModifier))
#endif
  {
    if(activeMdiChild())
    {
      BaseAdapterVisuQt *bavQt =
          dynamic_cast< BaseAdapterVisuQt * >(activeMdiChild());
      Adapter::Viewer *viewer = bavQt->getViewer();

#ifdef FEVV_USE_AIF
      if(dynamic_cast< SimpleViewer< FEVV::MeshAIF > * >(viewer))
        on_actionOpen_AIFMesh_ADD_SPACE_TIME(
            dynamic_cast< SimpleViewer< FEVV::MeshAIF > * >(viewer));
#else
      QMessageBox::information(
          this, "", QObject::tr("actionOpen_AIFMesh : FEVV_USE_AIF is OFF !"));
#endif
    }
  }
  else if(QApplication::keyboardModifiers().testFlag(
              Qt::ShiftModifier)) // QApplication::queryKeyboardModifiers()
    on_actionOpen_AIFMesh_SPACE_TIME();
#ifdef __APPLE__
  // on macOS, the ControlModifier value corresponds to the Command keys on the
  // keyboard, and the MetaModifier value corresponds to the Control keys
  else if(QApplication::keyboardModifiers().testFlag(
              /*Qt::MetaModifier*/ Qt::ControlModifier))
#else
  else if(QApplication::keyboardModifiers().testFlag(Qt::ControlModifier))
#endif
    on_actionOpen_AIFMesh_NORMAL(true); // empty mesh
  else
    on_actionOpen_AIFMesh_NORMAL();
}

template< typename HalfedgeGraph >
void
FEVV::SimpleWindow::writeHG(FEVV::SimpleViewer< HalfedgeGraph > *viewer)
{
  if(!viewer)
    return;

  QString validExtensions = "OBJ Files (*.obj);;OFF files (*.off);;COFF files "
                            "(*.coff);;PLY files (*.ply);;MSH files (*.msh)";
  QString validVtkExtensions =
      "VTK Files (*.vtk);;VTP files (*.vtp);;VTU files (*.vtu)";

  QString allExtensions = validExtensions
#ifdef FEVV_USE_VTK
                          + ";;" + validVtkExtensions
#endif
                          + "";

  std::vector< HalfedgeGraph * > meshes = viewer->getSelectedMeshes();
  std::vector< std::string > meshes_names = viewer->getSelectedMeshesNames();
  std::vector< FEVV::PMapsContainer * > properties_maps =
      viewer->getSelected_properties_maps();

  for(unsigned i = 0; i < meshes.size(); i++)
  {
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
#endif

    if(!fileName.isEmpty())
      FEVV::Filters::write_mesh(
          fileName.toStdString(), *(meshes[i]), *(properties_maps[i]));
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

#ifdef FEVV_USE_CGAL
      if(dynamic_cast< SimpleViewer< FEVV::MeshPolyhedron > * >(
             bav->getViewer()))
        writeHG(dynamic_cast< SimpleViewer< FEVV::MeshPolyhedron > * >(
            bav->getViewer()));
      else if(dynamic_cast< SimpleViewer< FEVV::MeshSurface > * >(
                  bav->getViewer()))
        writeHG(dynamic_cast< SimpleViewer< FEVV::MeshSurface > * >(
            bav->getViewer()));
      else if(dynamic_cast< SimpleViewer< FEVV::MeshLCC > * >(bav->getViewer()))
        writeHG(
            dynamic_cast< SimpleViewer< FEVV::MeshLCC > * >(bav->getViewer()));
#endif
#ifdef FEVV_USE_OPENMESH
      if(dynamic_cast< SimpleViewer< FEVV::MeshOpenMesh > * >(bav->getViewer()))
        writeHG(dynamic_cast< SimpleViewer< FEVV::MeshOpenMesh > * >(
            bav->getViewer()));
#endif
#ifdef FEVV_USE_AIF
      if(dynamic_cast< SimpleViewer< FEVV::MeshAIF > * >(bav->getViewer()))
        writeHG(
            dynamic_cast< SimpleViewer< FEVV::MeshAIF > * >(bav->getViewer()));
#endif
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
         "Open   -> <b>NORMAL</b>         mode : <b>no key</b><br>"
         "Open   -> <b>SPACE/TIME</b>     mode : <b>shift</b><br>"
         "Open   -> <b>EMPTY</b>          mode : <b>ctrl</b> *<br>"
         "Open   -> <b>ADD SPACE/TIME</b> mode : <b>shift+ctrl</b> *<br>"
         "<br>"
         "Viewer -> <b>SELECT</b>         mesh : <b>shift</b><br>"
         "Viewer -> <b>TRANSLATE</b>      mesh : <b>T</b> (Translation draggers must be shown, see toolbar)<br>"
         "Viewer -> <b>ROTATE</b>         mesh : <b>R</b> (Rotation draggers must be shown, see toolbar)<br>"
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

template< typename HalfedgeGraph >
void
FEVV::SimpleWindow::actionHG(FEVV::SimpleViewer< HalfedgeGraph > *viewer,
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

        //actionHG(viewer, 'D', '_'); // 19/03/19 - test for calling a re-'D'raw (because of 'hide/not hide' meshes in SPACE mode)
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
  else if(t == 'D') // 'D'raw
  {
    std::vector< HalfedgeGraph * > meshes = viewer->getMeshes();
    std::vector< std::string > meshes_names = viewer->getMeshesNames();
    std::vector< FEVV::PMapsContainer * > properties_maps =
        viewer->get_properties_maps();

    for(unsigned i = 0; i < meshes.size(); i++)
      draw_or_redraw_mesh(
          meshes[i], properties_maps[i], viewer, true, false, meshes_names[i]);
  }
}

inline void
FEVV::SimpleWindow::pre_actionHG(Adapter::Viewer *viewer, char t, char t2)
{
#ifdef FEVV_USE_CGAL
  if(dynamic_cast< SimpleViewer< FEVV::MeshPolyhedron > * >(viewer))
    actionHG(
        dynamic_cast< SimpleViewer< FEVV::MeshPolyhedron > * >(viewer), t, t2);
  else if(dynamic_cast< SimpleViewer< FEVV::MeshSurface > * >(viewer))
    actionHG(
        dynamic_cast< SimpleViewer< FEVV::MeshSurface > * >(viewer), t, t2);
  else if(dynamic_cast< SimpleViewer< FEVV::MeshLCC > * >(viewer))
    actionHG(dynamic_cast< SimpleViewer< FEVV::MeshLCC > * >(viewer), t, t2);
#endif
#ifdef FEVV_USE_OPENMESH
  if(dynamic_cast< SimpleViewer< FEVV::MeshOpenMesh > * >(viewer))
    actionHG(
        dynamic_cast< SimpleViewer< FEVV::MeshOpenMesh > * >(viewer), t, t2);
#endif
#ifdef FEVV_USE_AIF
  if(dynamic_cast< SimpleViewer< FEVV::MeshAIF > * >(viewer))
    actionHG(dynamic_cast< SimpleViewer< FEVV::MeshAIF > * >(viewer), t, t2);
#endif
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
    //bool SAVE_m_recreateOSGobj_if_redraw = bavQt->getViewer()->m_recreateOSGobj_if_redraw; bavQt->getViewer()->m_recreateOSGobj_if_redraw = true; // NEW
    pre_actionHG(bavQt->getViewer());
    //bavQt->getViewer()->m_recreateOSGobj_if_redraw = SAVE_m_recreateOSGobj_if_redraw; // NEW
    bavQt->getViewer()->m_space_time_changeColorMode = true;

    updateActiveChildTitle();
#endif
  }
}

template< typename HalfedgeGraph >
void
FEVV::SimpleWindow::centerHG(FEVV::SimpleViewer< HalfedgeGraph > *viewer)
{
  if(!viewer)
    return;

  std::vector< HalfedgeGraph * > meshes = viewer->getSelectedMeshes();

  for(unsigned i = 0; i < meshes.size(); i++)
  {
    viewer->centerMesh(meshes[i]);
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

#ifdef FEVV_USE_CGAL
      if(dynamic_cast< SimpleViewer< FEVV::MeshPolyhedron > * >(
             bav->getViewer()))
        centerHG(dynamic_cast< SimpleViewer< FEVV::MeshPolyhedron > * >(
            bav->getViewer()));
      else if(dynamic_cast< SimpleViewer< FEVV::MeshSurface > * >(
                  bav->getViewer()))
        centerHG(dynamic_cast< SimpleViewer< FEVV::MeshSurface > * >(
            bav->getViewer()));
      else if(dynamic_cast< SimpleViewer< FEVV::MeshLCC > * >(bav->getViewer()))
        centerHG(
            dynamic_cast< SimpleViewer< FEVV::MeshLCC > * >(bav->getViewer()));
#endif
#ifdef FEVV_USE_OPENMESH
      if(dynamic_cast< SimpleViewer< FEVV::MeshOpenMesh > * >(bav->getViewer()))
        centerHG(dynamic_cast< SimpleViewer< FEVV::MeshOpenMesh > * >(
            bav->getViewer()));
#endif
#ifdef FEVV_USE_AIF
      if(dynamic_cast< SimpleViewer< FEVV::MeshAIF > * >(bav->getViewer()))
        centerHG(
            dynamic_cast< SimpleViewer< FEVV::MeshAIF > * >(bav->getViewer()));
#endif
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

template< typename HalfedgeGraph >
void
FEVV::SimpleWindow::showSelectedHG(FEVV::SimpleViewer< HalfedgeGraph > *viewer)
{
  if(!viewer)
    return;
}

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
