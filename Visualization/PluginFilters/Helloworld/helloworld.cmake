# --> HelloworldPlugin : QtPlugin [BEGIN]
OPTION( BUILD_USE_GUI_HelloworldPlugin "BUILD HelloworldPlugin " ON )
if (BUILD_USE_GUI_HelloworldPlugin)
  set( Helloworld_Qt_Plugin_HEADER
       "${CMAKE_CURRENT_LIST_DIR}/helloworld_plugin.h"
       "${CMAKE_CURRENT_LIST_DIR}/Dialogs/helloworld_dialog.h"
       )
  set( Helloworld_Qt_Plugin_UI
       "${CMAKE_CURRENT_LIST_DIR}/Dialogs/helloworld_dialog.ui"
       )
  if (BUILD_USE_QT5)
    QT5_WRAP_CPP(Helloworld_Qt_Plugin_MOC_CPP ${Helloworld_Qt_Plugin_HEADER})
    QT5_WRAP_UI(Helloworld_Qt_Plugin_UI_CPP ${Helloworld_Qt_Plugin_UI})
    set(Helloworld_Qt_Plugin_SRC ${Helloworld_Qt_Plugin_SRC} ${Helloworld_Qt_Plugin_MOC_CPP} ${Helloworld_Qt_Plugin_UI_CPP})
  else(BUILD_USE_QT5)
    QT4_WRAP_CPP(Helloworld_Qt_Plugin_MOC_CPP ${Helloworld_Qt_Plugin_HEADER})
    QT4_WRAP_UI(Helloworld_Qt_Plugin_UI_CPP ${Helloworld_Qt_Plugin_UI})
    set(Helloworld_Qt_Plugin_SRC ${Helloworld_Qt_Plugin_SRC} ${Helloworld_Qt_Plugin_MOC_CPP} ${Helloworld_Qt_Plugin_UI_CPP})
  endif(BUILD_USE_QT5)

  add_library(
      HelloworldPlugin
      SHARED
      "${CMAKE_CURRENT_LIST_DIR}/helloworld_plugin.cpp"
      "${CMAKE_CURRENT_LIST_DIR}/Dialogs/helloworld_dialog.cpp"
      ${Helloworld_Qt_Plugin_SRC}
      ${osgQt_SRC} # from viewer
      )
  target_link_libraries (HelloworldPlugin ${Helloworld_Qt_Plugin_LIB}
    ${GUILIB_DEMO} # from viewer
    )
  add_dependencies(HelloworldPlugin mepp-gui)
endif (BUILD_USE_GUI_HelloworldPlugin)
# --> HelloworldPlugin : QtPlugin [END]

# remove the lines below to integrate into MEPP2-public
#TODO-elo why is this line needed ?
include_directories( "${CMAKE_CURRENT_BINARY_DIR}" )
