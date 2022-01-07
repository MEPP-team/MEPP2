# --> ProgressiveCompressionPlugin : QtPlugin [BEGIN]
if( NOT BUILD_USE_CGAL )
  return()
endif()
if( NOT BUILD_USE_DRACO )
  return()
endif()
OPTION( BUILD_USE_GUI_ProgressiveCompressionPlugin "BUILD ProgressiveCompressionPlugin " ON )
if (BUILD_USE_GUI_ProgressiveCompressionPlugin)
  set( ProgressiveCompression_Qt_Plugin_HEADER
       "${CMAKE_CURRENT_LIST_DIR}/progressivecompression_plugin.h"
       "${CMAKE_CURRENT_LIST_DIR}/Dialogs/progressivecompression_dialog.h"
       )
  set( ProgressiveCompression_Qt_Plugin_UI
       "${CMAKE_CURRENT_LIST_DIR}/Dialogs/progressivecompression_dialog.ui"
       )
  if (BUILD_USE_QT5)
    QT5_WRAP_CPP(ProgressiveCompression_Qt_Plugin_MOC_CPP ${ProgressiveCompression_Qt_Plugin_HEADER})
    QT5_WRAP_UI(ProgressiveCompression_Qt_Plugin_UI_CPP ${ProgressiveCompression_Qt_Plugin_UI})
    set(ProgressiveCompression_Qt_Plugin_SRC ${ProgressiveCompression_Qt_Plugin_SRC} ${ProgressiveCompression_Qt_Plugin_MOC_CPP} ${ProgressiveCompression_Qt_Plugin_UI_CPP})
  else(BUILD_USE_QT5)
    QT4_WRAP_CPP(ProgressiveCompression_Qt_Plugin_MOC_CPP ${ProgressiveCompression_Qt_Plugin_HEADER})
    QT4_WRAP_UI(ProgressiveCompression_Qt_Plugin_UI_CPP ${ProgressiveCompression_Qt_Plugin_UI})
    set(ProgressiveCompression_Qt_Plugin_SRC ${ProgressiveCompression_Qt_Plugin_SRC} ${ProgressiveCompression_Qt_Plugin_MOC_CPP} ${ProgressiveCompression_Qt_Plugin_UI_CPP})
  endif(BUILD_USE_QT5)
include_directories( ${draco_INCLUDE_DIRS} )
  add_library(
      ProgressiveCompressionPlugin
      SHARED
      "${CMAKE_CURRENT_LIST_DIR}/progressivecompression_plugin.cpp"
      "${CMAKE_CURRENT_LIST_DIR}/Dialogs/progressivecompression_dialog.cpp"
      ${ProgressiveCompression_Qt_Plugin_SRC}
      ${osgQt_SRC} # from viewer
      )
  target_link_libraries (ProgressiveCompressionPlugin ${ProgressiveCompression_Qt_Plugin_LIB}
    ${GUILIB_DEMO} # from viewer
    )
  add_dependencies(ProgressiveCompressionPlugin mepp-gui)
  target_link_libraries(ProgressiveCompressionPlugin ${draco_LIBRARIES})
endif (BUILD_USE_GUI_ProgressiveCompressionPlugin)
# --> ProgressiveCompressionPlugin : QtPlugin [END]

# remove the lines below to integrate into MEPP2-public
#TODO-elo why is this line needed ?
include_directories( "${CMAKE_CURRENT_BINARY_DIR}" )
