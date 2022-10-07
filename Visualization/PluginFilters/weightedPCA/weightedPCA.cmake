# The weightedPCA filter relies upon Flann,
# so the plugin should not appear when PCL is
# disabled
if( NOT BUILD_USE_PCL )
  return()
endif()

if( MSVC )
  set( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /IGNORE:4098" )
endif()

# --> WeightedPCAPlugin : QtPlugin [BEGIN]
OPTION( BUILD_USE_GUI_WeightedPCAPlugin "BUILD WeightedPCAPlugin " ON )
if (BUILD_USE_GUI_WeightedPCAPlugin)
  set( WeightedPCA_Qt_Plugin_HEADER
       "${CMAKE_CURRENT_LIST_DIR}/weightedPCA_plugin.h"
       "${CMAKE_CURRENT_LIST_DIR}/Dialogs/weightedPCA_dialog.h"
       )
  set( WeightedPCA_Qt_Plugin_UI
       "${CMAKE_CURRENT_LIST_DIR}/Dialogs/weightedPCA_dialog.ui"
       )
  if (BUILD_USE_QT6)
    QT6_WRAP_CPP(WeightedPCA_Qt_Plugin_MOC_CPP ${WeightedPCA_Qt_Plugin_HEADER})
    QT6_WRAP_UI(WeightedPCA_Qt_Plugin_UI_CPP ${WeightedPCA_Qt_Plugin_UI})
    set(WeightedPCA_Qt_Plugin_SRC ${WeightedPCA_Qt_Plugin_SRC} ${WeightedPCA_Qt_Plugin_MOC_CPP} ${WeightedPCA_Qt_Plugin_UI_CPP})
  elseif (BUILD_USE_QT5)
    QT5_WRAP_CPP(WeightedPCA_Qt_Plugin_MOC_CPP ${WeightedPCA_Qt_Plugin_HEADER})
    QT5_WRAP_UI(WeightedPCA_Qt_Plugin_UI_CPP ${WeightedPCA_Qt_Plugin_UI})
    set(WeightedPCA_Qt_Plugin_SRC ${WeightedPCA_Qt_Plugin_SRC} ${WeightedPCA_Qt_Plugin_MOC_CPP} ${WeightedPCA_Qt_Plugin_UI_CPP})
  else()
    QT4_WRAP_CPP(WeightedPCA_Qt_Plugin_MOC_CPP ${WeightedPCA_Qt_Plugin_HEADER})
    QT4_WRAP_UI(WeightedPCA_Qt_Plugin_UI_CPP ${WeightedPCA_Qt_Plugin_UI})
    set(WeightedPCA_Qt_Plugin_SRC ${WeightedPCA_Qt_Plugin_SRC} ${WeightedPCA_Qt_Plugin_MOC_CPP} ${WeightedPCA_Qt_Plugin_UI_CPP})
  endif()

  add_library(
      WeightedPCAPlugin
      SHARED
      "${CMAKE_CURRENT_LIST_DIR}/weightedPCA_plugin.cpp"
      "${CMAKE_CURRENT_LIST_DIR}/Dialogs/weightedPCA_dialog.cpp"
      ${WeightedPCA_Qt_Plugin_SRC}
      ${osgQt_SRC} ${MOC_FILES_osgQOpenGL} # from viewer
      )
  target_link_libraries (WeightedPCAPlugin ${WeightedPCA_Qt_Plugin_LIB}
    ${GUILIB_DEMO} # from viewer
    )
  add_dependencies(WeightedPCAPlugin mepp-gui)
endif (BUILD_USE_GUI_WeightedPCAPlugin)
# --> WeightedPCAPlugin : QtPlugin [END]

# remove the lines below to integrate into MEPP2-public
#TODO-elo why is this line needed ?
include_directories( "${CMAKE_CURRENT_BINARY_DIR}" )
