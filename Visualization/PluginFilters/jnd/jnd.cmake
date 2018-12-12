# --> JndPlugin : QtPlugin [BEGIN]
OPTION(BUILD_USE_GUI_JndPlugin "BUILD JndPlugin " ON)
if (BUILD_USE_GUI_JndPlugin)
  set(Jnd_Qt_Plugin_HEADER "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/jnd/JndPlugin.h" "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/jnd/Dialogs/DialogJnd1.h")
  set(Jnd_Qt_Plugin_UI "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/jnd/Dialogs/DialogJnd1.ui")
  if (BUILD_USE_QT5)
    QT5_WRAP_CPP(Jnd_Qt_Plugin_MOC_CPP ${Jnd_Qt_Plugin_HEADER})
    QT5_WRAP_UI(Jnd_Qt_Plugin_UI_CPP ${Jnd_Qt_Plugin_UI})
    set(Jnd_Qt_Plugin_SRC ${Jnd_Qt_Plugin_SRC} ${Jnd_Qt_Plugin_MOC_CPP} ${Jnd_Qt_Plugin_UI_CPP})
  else(BUILD_USE_QT5)
    QT4_WRAP_CPP(Jnd_Qt_Plugin_MOC_CPP ${Jnd_Qt_Plugin_HEADER})
    QT4_WRAP_UI(Jnd_Qt_Plugin_UI_CPP ${Jnd_Qt_Plugin_UI})
    set(Jnd_Qt_Plugin_SRC ${Jnd_Qt_Plugin_SRC} ${Jnd_Qt_Plugin_MOC_CPP} ${Jnd_Qt_Plugin_UI_CPP})
  endif(BUILD_USE_QT5)

  add_library(JndPlugin SHARED "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/jnd/JndPlugin.cpp" "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/jnd/Dialogs/DialogJnd1.cpp"
    ${Jnd_Qt_Plugin_SRC}
    ${osgQt_SRC} # from viewer
    )
  target_link_libraries (JndPlugin ${Jnd_Qt_Plugin_LIB}
    ${GUILIB_DEMO} # from viewer
    )
  add_dependencies(JndPlugin mepp-gui)
endif (BUILD_USE_GUI_JndPlugin)
# --> JndPlugin : QtPlugin [END]
