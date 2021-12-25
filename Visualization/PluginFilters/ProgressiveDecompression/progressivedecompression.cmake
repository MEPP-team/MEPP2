# --> ProgressiveCompressionPlugin : QtPlugin [BEGIN]
if( NOT BUILD_USE_CGAL )
  return()
endif()
if( NOT BUILD_USE_DRACO )
  return()
endif()
OPTION( BUILD_USE_GUI_ProgressiveDecompressionPlugin "BUILD ProgressiveDecompressionPlugin " ON )
if (BUILD_USE_GUI_ProgressiveDecompressionPlugin)
  set( ProgressiveDecompression_Qt_Plugin_HEADER
       "${CMAKE_CURRENT_LIST_DIR}/progressivedecompression_plugin.h"
       "${CMAKE_CURRENT_LIST_DIR}/Dialogs/progressivedecompression_dialog.h"
       )
  set( ProgressiveDecompression_Qt_Plugin_UI
       "${CMAKE_CURRENT_LIST_DIR}/Dialogs/progressivedecompression_dialog.ui"
       )
  if (BUILD_USE_QT5)
    QT5_WRAP_CPP(ProgressiveDecompression_Qt_Plugin_MOC_CPP ${ProgressiveDecompression_Qt_Plugin_HEADER})
    QT5_WRAP_UI(ProgressiveDecompression_Qt_Plugin_UI_CPP ${ProgressiveDecompression_Qt_Plugin_UI})
    set(ProgressiveDecompression_Qt_Plugin_SRC ${ProgressiveDecompression_Qt_Plugin_SRC} ${ProgressiveDecompression_Qt_Plugin_MOC_CPP} ${ProgressiveDecompression_Qt_Plugin_UI_CPP})
  else(BUILD_USE_QT5)
    QT4_WRAP_CPP(ProgressiveDecompression_Qt_Plugin_MOC_CPP ${ProgressiveDecompression_Qt_Plugin_HEADER})
    QT4_WRAP_UI(ProgressiveDecompression_Qt_Plugin_UI_CPP ${ProgressiveDecompression_Qt_Plugin_UI})
    set(ProgressiveDecompression_Qt_Plugin_SRC ${ProgressiveDecompression_Qt_Plugin_SRC} ${ProgressiveDecompression_Qt_Plugin_MOC_CPP} ${ProgressiveDecompression_Qt_Plugin_UI_CPP})
  endif(BUILD_USE_QT5)
include_directories( ${draco_INCLUDE_DIRS} )
  add_library(
      ProgressiveDecompressionPlugin
      SHARED
      "${CMAKE_CURRENT_LIST_DIR}/progressivedecompression_plugin.cpp"
      "${CMAKE_CURRENT_LIST_DIR}/Dialogs/progressivedecompression_dialog.cpp"
      ${ProgressiveDecompression_Qt_Plugin_SRC}
      ${osgQt_SRC} # from viewer
      )
  target_link_libraries (ProgressiveDecompressionPlugin ${ProgressiveDecompression_Qt_Plugin_LIB}
    ${GUILIB_DEMO} # from viewer
    )
  add_dependencies(ProgressiveDecompressionPlugin mepp-gui)
  target_link_libraries(ProgressiveDecompressionPlugin ${draco_LIBRARIES})
endif (BUILD_USE_GUI_ProgressiveDecompressionPlugin)
# --> ProgressiveDecompressionPlugin : QtPlugin [END]

# remove the lines below to integrate into MEPP2-public
#TODO-elo why is this line needed ?
include_directories( "${CMAKE_CURRENT_BINARY_DIR}" )
