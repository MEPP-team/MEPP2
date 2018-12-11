# --> MSDM2Plugin : QtPlugin [BEGIN]
OPTION(BUILD_USE_GUI_MSDM2Plugin "BUILD MSDM2Plugin " ON)
if (BUILD_USE_GUI_MSDM2Plugin AND BUILD_USE_CGAL)
  set(MSDM2_Qt_Plugin_HEADER "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/msdm2/MSDM2Plugin.h" "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/msdm2/Dialogs/DialogMSDM21.h")
  set(MSDM2_Qt_Plugin_UI "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/msdm2/Dialogs/DialogMSDM21.ui")
  if (BUILD_USE_QT5)
    QT5_WRAP_CPP(MSDM2_Qt_Plugin_MOC_CPP ${MSDM2_Qt_Plugin_HEADER})
    QT5_WRAP_UI(MSDM2_Qt_Plugin_UI_CPP ${MSDM2_Qt_Plugin_UI})
    set(MSDM2_Qt_Plugin_SRC ${MSDM2_Qt_Plugin_SRC} ${MSDM2_Qt_Plugin_MOC_CPP} ${MSDM2_Qt_Plugin_UI_CPP})
  else(BUILD_USE_QT5)
    QT4_WRAP_CPP(MSDM2_Qt_Plugin_MOC_CPP ${MSDM2_Qt_Plugin_HEADER})
    QT4_WRAP_UI(MSDM2_Qt_Plugin_UI_CPP ${MSDM2_Qt_Plugin_UI})
    set(MSDM2_Qt_Plugin_SRC ${MSDM2_Qt_Plugin_SRC} ${MSDM2_Qt_Plugin_MOC_CPP} ${MSDM2_Qt_Plugin_UI_CPP})
  endif(BUILD_USE_QT5)

  add_library(MSDM2Plugin SHARED
    "${PROJECT_SOURCE_DIR}/FEVV/Filters/CGAL/Surface_mesh/msdm2.h"
    "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/msdm2/MSDM2Plugin.cpp"
    "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/msdm2/Dialogs/DialogMSDM21.cpp"
    ${MSDM2_Qt_Plugin_SRC}
    ${osgQt_SRC} # from viewer
    )
  target_link_libraries (MSDM2Plugin ${MSDM2_Qt_Plugin_LIB}
    ${GUILIB_DEMO} # from viewer
    )
  add_dependencies(MSDM2Plugin mepp-gui)
endif ()
# --> MSDM2Plugin : QtPlugin [END]
