# --> CopyGraphPlugin : QtPlugin [BEGIN]
OPTION( BUILD_USE_GUI_CopyGraphPlugin "BUILD CopyGraphPlugin " ON )
if (BUILD_USE_GUI_CopyGraphPlugin)
  set( CopyGraph_Qt_Plugin_HEADER
       "${CMAKE_CURRENT_LIST_DIR}/copy_graph_plugin.h"
       "${CMAKE_CURRENT_LIST_DIR}/Dialogs/copy_graph_dialog.h"
       )
  set( CopyGraph_Qt_Plugin_UI
       "${CMAKE_CURRENT_LIST_DIR}/Dialogs/copy_graph_dialog.ui"
       )
  if (BUILD_USE_QT5)
    QT5_WRAP_CPP(CopyGraph_Qt_Plugin_MOC_CPP ${CopyGraph_Qt_Plugin_HEADER})
    QT5_WRAP_UI(CopyGraph_Qt_Plugin_UI_CPP ${CopyGraph_Qt_Plugin_UI})
    set(CopyGraph_Qt_Plugin_SRC ${CopyGraph_Qt_Plugin_SRC} ${CopyGraph_Qt_Plugin_MOC_CPP} ${CopyGraph_Qt_Plugin_UI_CPP})
  else(BUILD_USE_QT5)
    QT4_WRAP_CPP(CopyGraph_Qt_Plugin_MOC_CPP ${CopyGraph_Qt_Plugin_HEADER})
    QT4_WRAP_UI(CopyGraph_Qt_Plugin_UI_CPP ${CopyGraph_Qt_Plugin_UI})
    set(CopyGraph_Qt_Plugin_SRC ${CopyGraph_Qt_Plugin_SRC} ${CopyGraph_Qt_Plugin_MOC_CPP} ${CopyGraph_Qt_Plugin_UI_CPP})
  endif(BUILD_USE_QT5)

  add_library(
      CopyGraphPlugin
      SHARED
      "${CMAKE_CURRENT_LIST_DIR}/copy_graph_plugin.cpp"
      "${CMAKE_CURRENT_LIST_DIR}/Dialogs/copy_graph_dialog.cpp"
      ${CopyGraph_Qt_Plugin_SRC}
      ${osgQt_SRC} ${MOC_FILES_osgQOpenGL} # from viewer
      )
  target_link_libraries (CopyGraphPlugin ${CopyGraph_Qt_Plugin_LIB}
    ${GUILIB_DEMO} # from viewer
    )
  add_dependencies(CopyGraphPlugin mepp-gui)
endif (BUILD_USE_GUI_CopyGraphPlugin)
# --> CopyGraphPlugin : QtPlugin [END]

# remove the lines below to integrate into MEPP2-public
#TODO-elo why is this line needed ?
include_directories( "${CMAKE_CURRENT_BINARY_DIR}" )
