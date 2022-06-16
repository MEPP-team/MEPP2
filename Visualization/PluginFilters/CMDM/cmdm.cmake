# the Boolean Operations filter relies upon CGAL,
# so the plugin should not appear when CGAL is
# disabled
if( NOT BUILD_USE_CGAL )
  return()
endif()

# --> CMDMPlugin : QtPlugin [BEGIN]
OPTION( BUILD_USE_GUI_CMDMPlugin "BUILD CMDMPlugin " ON )
if (BUILD_USE_GUI_CMDMPlugin)
  set( CMDM_Qt_Plugin_HEADER
       "${CMAKE_CURRENT_LIST_DIR}/cmdm_plugin.h"
       "${CMAKE_CURRENT_LIST_DIR}/Dialogs/cmdm_dialog.h"
       )
  set( CMDM_Qt_Plugin_UI
       "${CMAKE_CURRENT_LIST_DIR}/Dialogs/cmdm_dialog.ui"
       )
  if (BUILD_USE_QT5)
    QT5_WRAP_CPP(CMDM_Qt_Plugin_MOC_CPP ${CMDM_Qt_Plugin_HEADER})
    QT5_WRAP_UI(CMDM_Qt_Plugin_UI_CPP ${CMDM_Qt_Plugin_UI})
    set(CMDM_Qt_Plugin_SRC ${CMDM_Qt_Plugin_SRC} ${CMDM_Qt_Plugin_MOC_CPP} ${CMDM_Qt_Plugin_UI_CPP})
  else(BUILD_USE_QT5)
    QT4_WRAP_CPP(CMDM_Qt_Plugin_MOC_CPP ${CMDM_Qt_Plugin_HEADER})
    QT4_WRAP_UI(CMDM_Qt_Plugin_UI_CPP ${CMDM_Qt_Plugin_UI})
    set(CMDM_Qt_Plugin_SRC ${CMDM_Qt_Plugin_SRC} ${CMDM_Qt_Plugin_MOC_CPP} ${CMDM_Qt_Plugin_UI_CPP})
  endif(BUILD_USE_QT5)

  add_library(
      CMDMPlugin
      SHARED
      "${CMAKE_CURRENT_LIST_DIR}/cmdm_plugin.cpp"
      "${CMAKE_CURRENT_LIST_DIR}/Dialogs/cmdm_dialog.cpp"
      ${CMDM_Qt_Plugin_SRC}
      ${osgQt_SRC} ${MOC_FILES_osgQOpenGL} # from viewer
      )
  target_link_libraries (CMDMPlugin ${CMDM_Qt_Plugin_LIB}
    ${GUILIB_DEMO} # from viewer
    )
  add_dependencies(CMDMPlugin mepp-gui)
endif (BUILD_USE_GUI_CMDMPlugin)
# --> CMDMPlugin : QtPlugin [END]

# remove the lines below to integrate into MEPP2-public
#TODO-elo why is this line needed ?
include_directories( "${CMAKE_CURRENT_BINARY_DIR}" )
