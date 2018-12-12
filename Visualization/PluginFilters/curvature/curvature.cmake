# --> CurvaturePlugin : QtPlugin [BEGIN]
OPTION(BUILD_USE_GUI_CurvaturePlugin "BUILD CurvaturePlugin " ON)
if (BUILD_USE_GUI_CurvaturePlugin)
  set(Curvature_Qt_Plugin_HEADER "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/curvature/CurvaturePlugin.h" "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/curvature/Dialogs/DialogCurvature1.h")
  set(Curvature_Qt_Plugin_UI "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/curvature/Dialogs/DialogCurvature1.ui")
  if (BUILD_USE_QT5)
    QT5_WRAP_CPP(Curvature_Qt_Plugin_MOC_CPP ${Curvature_Qt_Plugin_HEADER})
    QT5_WRAP_UI(Curvature_Qt_Plugin_UI_CPP ${Curvature_Qt_Plugin_UI})
    set(Curvature_Qt_Plugin_SRC ${Curvature_Qt_Plugin_SRC} ${Curvature_Qt_Plugin_MOC_CPP} ${Curvature_Qt_Plugin_UI_CPP})
  else(BUILD_USE_QT5)
    QT4_WRAP_CPP(Curvature_Qt_Plugin_MOC_CPP ${Curvature_Qt_Plugin_HEADER})
    QT4_WRAP_UI(Curvature_Qt_Plugin_UI_CPP ${Curvature_Qt_Plugin_UI})
    set(Curvature_Qt_Plugin_SRC ${Curvature_Qt_Plugin_SRC} ${Curvature_Qt_Plugin_MOC_CPP} ${Curvature_Qt_Plugin_UI_CPP})
  endif(BUILD_USE_QT5)

  add_library(CurvaturePlugin SHARED "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/curvature/CurvaturePlugin.cpp" "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/curvature/Dialogs/DialogCurvature1.cpp"
    ${Curvature_Qt_Plugin_SRC}
    ${osgQt_SRC} # from viewer
    )
  target_link_libraries (CurvaturePlugin ${Curvature_Qt_Plugin_LIB}
    ${GUILIB_DEMO} # from viewer
    )
  add_dependencies(CurvaturePlugin mepp-gui)
endif (BUILD_USE_GUI_CurvaturePlugin)
# --> CurvaturePlugin : QtPlugin [END]
