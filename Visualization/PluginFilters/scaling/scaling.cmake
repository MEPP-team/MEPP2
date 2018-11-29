# --> ScalingPlugin : QtPlugin [BEGIN]
OPTION(BUILD_USE_GUI_ScalingPlugin "BUILD ScalingPlugin " ON)
if (BUILD_USE_GUI_ScalingPlugin)
  set(Scaling_Qt_Plugin_HEADER "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/scaling/ScalingPlugin.h" "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/scaling/Dialogs/DialogScaling1.h")
  set(Scaling_Qt_Plugin_UI "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/scaling/Dialogs/DialogScaling1.ui")
  if (BUILD_USE_QT5)
    QT5_WRAP_CPP(Scaling_Qt_Plugin_MOC_CPP ${Scaling_Qt_Plugin_HEADER})
    QT5_WRAP_UI(Scaling_Qt_Plugin_UI_CPP ${Scaling_Qt_Plugin_UI})
    set(Scaling_Qt_Plugin_SRC ${Scaling_Qt_Plugin_SRC} ${Scaling_Qt_Plugin_MOC_CPP} ${Scaling_Qt_Plugin_UI_CPP})
  else(BUILD_USE_QT5)
    QT4_WRAP_CPP(Scaling_Qt_Plugin_MOC_CPP ${Scaling_Qt_Plugin_HEADER})
    QT4_WRAP_UI(Scaling_Qt_Plugin_UI_CPP ${Scaling_Qt_Plugin_UI})
    set(Scaling_Qt_Plugin_SRC ${Scaling_Qt_Plugin_SRC} ${Scaling_Qt_Plugin_MOC_CPP} ${Scaling_Qt_Plugin_UI_CPP})
  endif(BUILD_USE_QT5)

  add_library(ScalingPlugin SHARED "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/scaling/ScalingPlugin.cpp" "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/scaling/Dialogs/DialogScaling1.cpp"
    ${Scaling_Qt_Plugin_SRC}
    ${osgQt_SRC} # from viewer
    )
  target_link_libraries (ScalingPlugin ${Scaling_Qt_Plugin_LIB}
    ${GUILIB_DEMO} # from viewer
    )
  add_dependencies(ScalingPlugin mepp-gui)
endif (BUILD_USE_GUI_ScalingPlugin)
# --> ScalingPlugin : QtPlugin [END]
