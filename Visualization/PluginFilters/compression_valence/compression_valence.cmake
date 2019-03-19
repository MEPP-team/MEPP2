# --> CompressionValencePlugin : QtPlugin [BEGIN]
OPTION(BUILD_USE_GUI_CompressionValencePlugin "BUILD CompressionValencePlugin " ON)
if (BUILD_USE_GUI_CompressionValencePlugin)
  set(CompressionValence_Qt_Plugin_HEADER "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/compression_valence/CompressionValencePlugin.h" "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/compression_valence/Dialogs/DialogCompressionValence1.h")
  set(CompressionValence_Qt_Plugin_UI "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/compression_valence/Dialogs/DialogCompressionValence1.ui")
  if (BUILD_USE_QT5)
    QT5_WRAP_CPP(CompressionValence_Qt_Plugin_MOC_CPP ${CompressionValence_Qt_Plugin_HEADER})
    QT5_WRAP_UI(CompressionValence_Qt_Plugin_UI_CPP ${CompressionValence_Qt_Plugin_UI})
    set(CompressionValence_Qt_Plugin_SRC ${CompressionValence_Qt_Plugin_SRC} ${CompressionValence_Qt_Plugin_MOC_CPP} ${CompressionValence_Qt_Plugin_UI_CPP})
  else(BUILD_USE_QT5)
    QT4_WRAP_CPP(CompressionValence_Qt_Plugin_MOC_CPP ${CompressionValence_Qt_Plugin_HEADER})
    QT4_WRAP_UI(CompressionValence_Qt_Plugin_UI_CPP ${CompressionValence_Qt_Plugin_UI})
    set(CompressionValence_Qt_Plugin_SRC ${CompressionValence_Qt_Plugin_SRC} ${CompressionValence_Qt_Plugin_MOC_CPP} ${CompressionValence_Qt_Plugin_UI_CPP})
  endif(BUILD_USE_QT5)

  add_library(CompressionValencePlugin SHARED
              "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/compression_valence/CompressionValencePlugin.cpp"
              "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/compression_valence/Dialogs/DialogCompressionValence1.cpp"
              ${CompressionValence_Qt_Plugin_SRC}
              ${osgQt_SRC} # from viewer
              )
  target_link_libraries (CompressionValencePlugin ${CompressionValence_Qt_Plugin_LIB}
    ${GUILIB_DEMO} # from viewer
    )
  add_dependencies(CompressionValencePlugin mepp-gui)
endif (BUILD_USE_GUI_CompressionValencePlugin)
# --> CompressionValencePlugin : QtPlugin [END]
