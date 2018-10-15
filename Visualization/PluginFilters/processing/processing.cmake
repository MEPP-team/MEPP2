# --> ProcessingPlugin : QtPlugin [BEGIN]
OPTION(BUILD_USE_GUI_ProcessingPlugin "BUILD ProcessingPlugin " ON)
if (BUILD_USE_GUI_ProcessingPlugin)
  set(Processing_Qt_Plugin_HEADER "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/processing/ProcessingPlugin.h" "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/processing/Dialogs/DialogProcessing1.h")
  set(Processing_Qt_Plugin_UI "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/processing/Dialogs/DialogProcessing1.ui")
  if (BUILD_USE_QT5)
    QT5_WRAP_CPP(Processing_Qt_Plugin_MOC_CPP ${Processing_Qt_Plugin_HEADER})
    QT5_WRAP_UI(Processing_Qt_Plugin_UI_CPP ${Processing_Qt_Plugin_UI})
    set(Processing_Qt_Plugin_SRC ${Processing_Qt_Plugin_SRC} ${Processing_Qt_Plugin_MOC_CPP} ${Processing_Qt_Plugin_UI_CPP})
  else(BUILD_USE_QT5)
    QT4_WRAP_CPP(Processing_Qt_Plugin_MOC_CPP ${Processing_Qt_Plugin_HEADER})
    QT4_WRAP_UI(Processing_Qt_Plugin_UI_CPP ${Processing_Qt_Plugin_UI})
    set(Processing_Qt_Plugin_SRC ${Processing_Qt_Plugin_SRC} ${Processing_Qt_Plugin_MOC_CPP} ${Processing_Qt_Plugin_UI_CPP})
  endif(BUILD_USE_QT5)

  add_library(ProcessingPlugin SHARED "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/processing/ProcessingPlugin.cpp" "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/processing/Dialogs/DialogProcessing1.cpp"
    ${Processing_Qt_Plugin_SRC}
    ${osgQt_SRC} # from viewer
    ${GUISRC_PLUGIN}
    )
  target_link_libraries (ProcessingPlugin ${Processing_Qt_Plugin_LIB}
    ${GUILIB_DEMO} # from viewer
    )
  add_dependencies(ProcessingPlugin mepp)
endif (BUILD_USE_GUI_ProcessingPlugin)
# --> ProcessingPlugin : QtPlugin [END]
