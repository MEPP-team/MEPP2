# --> TextureImageDemoPlugin : QtPlugin [BEGIN]
OPTION( BUILD_USE_GUI_TextureImageDemoPlugin "BUILD TextureImageDemoPlugin " ON )
if (BUILD_USE_GUI_TextureImageDemoPlugin)
  set( TextureImageDemo_Qt_Plugin_HEADER
       "${CMAKE_CURRENT_LIST_DIR}/texture_image_demo_plugin.h"
       "${CMAKE_CURRENT_LIST_DIR}/Dialogs/texture_image_demo_dialog.h"
       )
  set( TextureImageDemo_Qt_Plugin_UI
       "${CMAKE_CURRENT_LIST_DIR}/Dialogs/texture_image_demo_dialog.ui"
       )
  if (BUILD_USE_QT5)
    QT5_WRAP_CPP(TextureImageDemo_Qt_Plugin_MOC_CPP ${TextureImageDemo_Qt_Plugin_HEADER})
    QT5_WRAP_UI(TextureImageDemo_Qt_Plugin_UI_CPP ${TextureImageDemo_Qt_Plugin_UI})
    set(TextureImageDemo_Qt_Plugin_SRC ${TextureImageDemo_Qt_Plugin_SRC} ${TextureImageDemo_Qt_Plugin_MOC_CPP} ${TextureImageDemo_Qt_Plugin_UI_CPP})
  else(BUILD_USE_QT5)
    QT4_WRAP_CPP(TextureImageDemo_Qt_Plugin_MOC_CPP ${TextureImageDemo_Qt_Plugin_HEADER})
    QT4_WRAP_UI(TextureImageDemo_Qt_Plugin_UI_CPP ${TextureImageDemo_Qt_Plugin_UI})
    set(TextureImageDemo_Qt_Plugin_SRC ${TextureImageDemo_Qt_Plugin_SRC} ${TextureImageDemo_Qt_Plugin_MOC_CPP} ${TextureImageDemo_Qt_Plugin_UI_CPP})
  endif(BUILD_USE_QT5)

  add_library(
      TextureImageDemoPlugin
      SHARED
      "${CMAKE_CURRENT_LIST_DIR}/texture_image_demo_plugin.cpp"
      "${CMAKE_CURRENT_LIST_DIR}/Dialogs/texture_image_demo_dialog.cpp"
      ${TextureImageDemo_Qt_Plugin_SRC}
      ${osgQt_SRC} ${MOC_FILES_osgQOpenGL} # from viewer
      )
  target_link_libraries (TextureImageDemoPlugin ${TextureImageDemo_Qt_Plugin_LIB}
    ${GUILIB_DEMO} # from viewer
    )
  add_dependencies(TextureImageDemoPlugin mepp-gui)
endif (BUILD_USE_GUI_TextureImageDemoPlugin)
# --> TextureImageDemoPlugin : QtPlugin [END]

# remove the lines below to integrate into MEPP2-public
#TODO-elo why is this line needed ?
include_directories( "${CMAKE_CURRENT_BINARY_DIR}" )
