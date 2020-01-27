# --> PointCloudCurvaturePlugin : QtPlugin [BEGIN]
OPTION( BUILD_USE_GUI_PointCloudCurvaturePlugin "BUILD PointCloudCurvaturePlugin " ON )
if (BUILD_USE_GUI_PointCloudCurvaturePlugin)
  set( PointCloudCurvature_Qt_Plugin_HEADER
       "${CMAKE_CURRENT_LIST_DIR}/point_cloud_curvature_plugin.h"
       "${CMAKE_CURRENT_LIST_DIR}/Dialogs/point_cloud_curvature_dialog.h"
       )
  set( PointCloudCurvature_Qt_Plugin_UI
       "${CMAKE_CURRENT_LIST_DIR}/Dialogs/point_cloud_curvature_dialog.ui"
       )
  if (BUILD_USE_QT5)
    QT5_WRAP_CPP(PointCloudCurvature_Qt_Plugin_MOC_CPP ${PointCloudCurvature_Qt_Plugin_HEADER})
    QT5_WRAP_UI(PointCloudCurvature_Qt_Plugin_UI_CPP ${PointCloudCurvature_Qt_Plugin_UI})
    set(PointCloudCurvature_Qt_Plugin_SRC ${PointCloudCurvature_Qt_Plugin_SRC} ${PointCloudCurvature_Qt_Plugin_MOC_CPP} ${PointCloudCurvature_Qt_Plugin_UI_CPP})
  else(BUILD_USE_QT5)
    QT4_WRAP_CPP(PointCloudCurvature_Qt_Plugin_MOC_CPP ${PointCloudCurvature_Qt_Plugin_HEADER})
    QT4_WRAP_UI(PointCloudCurvature_Qt_Plugin_UI_CPP ${PointCloudCurvature_Qt_Plugin_UI})
    set(PointCloudCurvature_Qt_Plugin_SRC ${PointCloudCurvature_Qt_Plugin_SRC} ${PointCloudCurvature_Qt_Plugin_MOC_CPP} ${PointCloudCurvature_Qt_Plugin_UI_CPP})
  endif(BUILD_USE_QT5)

  add_library(
      PointCloudCurvaturePlugin
      SHARED
      "${CMAKE_CURRENT_LIST_DIR}/point_cloud_curvature_plugin.cpp"
      "${CMAKE_CURRENT_LIST_DIR}/Dialogs/point_cloud_curvature_dialog.cpp"
      ${PointCloudCurvature_Qt_Plugin_SRC}
      ${osgQt_SRC} # from viewer
      )
  target_link_libraries (PointCloudCurvaturePlugin ${PointCloudCurvature_Qt_Plugin_LIB}
    ${GUILIB_DEMO} # from viewer
    )
  add_dependencies(PointCloudCurvaturePlugin mepp-gui)
endif (BUILD_USE_GUI_PointCloudCurvaturePlugin)
# --> PointCloudCurvaturePlugin : QtPlugin [END]

# remove the lines below to integrate into MEPP2-public
#TODO-elo why is this line needed ?
include_directories( "${CMAKE_CURRENT_BINARY_DIR}" )
