# --> DecompressionValencePlugin : QtPlugin [BEGIN]
OPTION(BUILD_USE_GUI_DecompressionValencePlugin "BUILD DecompressionValencePlugin " ON)
if (BUILD_USE_GUI_DecompressionValencePlugin)
  set(DecompressionValence_Qt_Plugin_HEADER "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/decompression_valence/DecompressionValencePlugin.h" "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/decompression_valence/Dialogs/DialogDecompressionValence1.h")
  set(DecompressionValence_Qt_Plugin_UI "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/decompression_valence/Dialogs/DialogDecompressionValence1.ui")
  if (BUILD_USE_QT5)
    QT5_WRAP_CPP(DecompressionValence_Qt_Plugin_MOC_CPP ${DecompressionValence_Qt_Plugin_HEADER})
    QT5_WRAP_UI(DecompressionValence_Qt_Plugin_UI_CPP ${DecompressionValence_Qt_Plugin_UI})
    set(DecompressionValence_Qt_Plugin_SRC ${DecompressionValence_Qt_Plugin_SRC} ${DecompressionValence_Qt_Plugin_MOC_CPP} ${DecompressionValence_Qt_Plugin_UI_CPP})
  else(BUILD_USE_QT5)
    QT4_WRAP_CPP(DecompressionValence_Qt_Plugin_MOC_CPP ${DecompressionValence_Qt_Plugin_HEADER})
    QT4_WRAP_UI(DecompressionValence_Qt_Plugin_UI_CPP ${DecompressionValence_Qt_Plugin_UI})
    set(DecompressionValence_Qt_Plugin_SRC ${DecompressionValence_Qt_Plugin_SRC} ${DecompressionValence_Qt_Plugin_MOC_CPP} ${DecompressionValence_Qt_Plugin_UI_CPP})
  endif(BUILD_USE_QT5)

  add_library(DecompressionValencePlugin SHARED
              "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/decompression_valence/DecompressionValencePlugin.cpp"
              "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/decompression_valence/Dialogs/DialogDecompressionValence1.cpp"
              ${DecompressionValence_Qt_Plugin_SRC}
              ${osgQt_SRC} # from viewer
              )
  target_link_libraries (DecompressionValencePlugin ${DecompressionValence_Qt_Plugin_LIB}
    ${GUILIB_DEMO} # from viewer
    )
  add_dependencies(DecompressionValencePlugin mepp-gui)
  if( MSVC )
    set_target_properties(DecompressionValencePlugin PROPERTIES LINK_FLAGS "/FORCE:MULTIPLE /IGNORE:4006 /IGNORE:4088")
  endif( MSVC )
endif (BUILD_USE_GUI_DecompressionValencePlugin)
# --> DecompressionValencePlugin : QtPlugin [END]
