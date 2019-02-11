# the Boolean Operations filter relies upon CGAL,
# so the plugin should not appear when CGAL is
# disabled
if( NOT BUILD_USE_CGAL )
  return()
endif()

# --> BooleanOperationsPlugin : QtPlugin [BEGIN]
OPTION(BUILD_USE_GUI_BooleanOperationsPlugin "BUILD BooleanOperationsPlugin " ON)
if (BUILD_USE_GUI_BooleanOperationsPlugin)
  set(BooleanOperations_Qt_Plugin_HEADER "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/boolean_operations/BooleanOperationsPlugin.h" "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/boolean_operations/Dialogs/DialogBooleanOperations1.h")
  set(BooleanOperations_Qt_Plugin_UI "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/boolean_operations/Dialogs/DialogBooleanOperations1.ui")
  if (BUILD_USE_QT5)
    QT5_WRAP_CPP(BooleanOperations_Qt_Plugin_MOC_CPP ${BooleanOperations_Qt_Plugin_HEADER})
    QT5_WRAP_UI(BooleanOperations_Qt_Plugin_UI_CPP ${BooleanOperations_Qt_Plugin_UI})
    set(BooleanOperations_Qt_Plugin_SRC ${BooleanOperations_Qt_Plugin_SRC} ${BooleanOperations_Qt_Plugin_MOC_CPP} ${BooleanOperations_Qt_Plugin_UI_CPP})
  else(BUILD_USE_QT5)
    QT4_WRAP_CPP(BooleanOperations_Qt_Plugin_MOC_CPP ${BooleanOperations_Qt_Plugin_HEADER})
    QT4_WRAP_UI(BooleanOperations_Qt_Plugin_UI_CPP ${BooleanOperations_Qt_Plugin_UI})
    set(BooleanOperations_Qt_Plugin_SRC ${BooleanOperations_Qt_Plugin_SRC} ${BooleanOperations_Qt_Plugin_MOC_CPP} ${BooleanOperations_Qt_Plugin_UI_CPP})
  endif(BUILD_USE_QT5)

  add_library(BooleanOperationsPlugin SHARED "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/boolean_operations/BooleanOperationsPlugin.cpp" "${PROJECT_SOURCE_DIR}/Visualization/PluginFilters/boolean_operations/Dialogs/DialogBooleanOperations1.cpp"
    ${BooleanOperations_Qt_Plugin_SRC}
    ${osgQt_SRC} # from viewer
    )
  target_link_libraries (BooleanOperationsPlugin ${BooleanOperations_Qt_Plugin_LIB}
    ${GUILIB_DEMO} # from viewer
    )
  add_dependencies(BooleanOperationsPlugin mepp-gui)
  if( MSVC )
    set_target_properties(BooleanOperationsPlugin PROPERTIES LINK_FLAGS "/FORCE:MULTIPLE /IGNORE:4006 /IGNORE:4088")
  endif( MSVC )
endif (BUILD_USE_GUI_BooleanOperationsPlugin)
# --> BooleanOperationsPlugin : QtPlugin [END]
