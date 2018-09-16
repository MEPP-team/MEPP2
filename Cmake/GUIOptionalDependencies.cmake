# -----------------------------------------------------------------------------
# Check Optional Dependencies
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Global options
# -----------------------------------------------------------------------------
message(STATUS "--------------------------------------------------------------")
message(STATUS "MEPP2 GUI optional configuration:")
message(STATUS "   (to change these values, use ccmake, a graphical")
message(STATUS "   cmake frontend, or define cmake commandline variables")
message(STATUS "   -e.g. '-DBUILD_USE_GUI:string=true'-, cf documentation)")
message(STATUS "")

option( BUILD_USE_QT5 "Using Qt5." OFF )

if( BUILD_USE_QT5 )
  set( LIST_OPTION ${LIST_OPTION} [QT5]\ )
  message(STATUS "   BUILD_USE_QT5  true   (Using of Qt5 instead of Qt4)")
else()
  message(STATUS "   BUILD_USE_QT5  false  (Using of Qt5 instead of Qt4)")
endif()

if( BUILD_USE_GUI )
  set( LIST_OPTION ${LIST_OPTION} [UI]\ )
  message(STATUS "   BUILD_USE_GUI  true   (Build user interface)")
else()
  message(STATUS "   BUILD_USE_GUI  false  (Don't build User interface (cli only))")
endif()

# -----------------------------------------------------------------------------
# Look for Qt (needed by openscenegraph visualization).
# -----------------------------------------------------------------------------
set( QT4_FOUND_MEPP 0 )
set( QT5_FOUND_MEPP 0 )
if( BUILD_USE_GUI )
  if( BUILD_USE_QT5 ) # sample for help : QT5_DIR=/usr/local/Cellar/qt/5.9.2 or Qt5_DIR=/usr/local/Cellar/qt/5.9.2/lib/cmake/Qt5 (caution, here 't' not 'T' !)
    # QT5 Handling
    if(DEFINED ENV{QT5_DIR})
      set( QT5_DIR $ENV{QT5_DIR} )
    endif()
    set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${QT5_DIR})

    # Instruct CMake to NOT run moc automatically when needed
    set(CMAKE_AUTOMOC OFF)

    find_package(Qt5 COMPONENTS Widgets OpenGL Xml REQUIRED)

    if( Qt5Widgets_FOUND AND Qt5OpenGL_FOUND AND Qt5Xml_FOUND )
      set(QT5_FOUND_MEPP 1)
      message( STATUS "Qt5 (Widgets, OpenGL and Xml modules) found (needed by OpenSceneGraph compiled with Qt5)." )

      add_definitions("-DFEVV_USE_QT5")

      set(MEPP_GUI_LIB ${MEPP_GUI_LIB}
        ${Qt5Widgets_LIBRARIES}
        ${Qt5OpenGL_LIBRARIES}
        ${Qt5Xml_LIBRARIES})
      set(MEPP_GUI_INC ${MEPP_GUI_INC}
        ${Qt5Widgets_INCLUDES_DIRS}
        ${Qt5OpenGL_INCLUDES_DIR}
        ${Qt5Xml_INCLUDES_DIR})

    else()
      message( STATUS "One of Qt5's modules was not found (needed by OpenSceneGraph)." )
    endif()
  else()
    # QT4 Handling
    if( QTDIR )
      set( ENV{QTDIR} ${QTDIR} )
    endif()
    set( QT_USE_QTMAIN    TRUE )
    set( QT_USE_QTXML     TRUE )
    set( QT_USE_QTOPENGL  TRUE )

    find_package(Qt4 COMPONENTS QtCore QtGUI QtXml QtOpenGL REQUIRED)

    if( QT4_FOUND )
      set(QT4_FOUND_MEPP 1)
      message(STATUS  "Qt4 found (needed by OpenSceneGraph).")

      add_definitions("-DFEVV_USE_QT4")

      set(QT_USE_QTXML 1)

      include(${QT_USE_FILE})

      set( MEPP_GUI_LIB ${MEPP_GUI_LIB} ${QT_LIBRARIES} )
      set( MEPP_GUI_INC ${MEPP_GUI_INC} ${QT_INCLUDE_DIR} )
    else()
      message(FATAL_ERROR "Qt4 not found (needed by OpenSceneGraph). Check the cmake variables associated to this package or disable it.")
    endif()
  endif()
endif()

# -----------------------------------------------------------------------------
# Look for OpenSceneGraph for 3D display.
# (They are not compulsory).
# -----------------------------------------------------------------------------
set( OPENSCENEGRAPH_FOUND_MEPP 0 )

if( BUILD_USE_GUI )

  # Find OpenSceneGraph
  find_package(OpenSceneGraph)
  if( OPENSCENEGRAPH_FOUND )
    set(MEPP_GUI_INC ${MEPP_GUI_INC} ${OPENSCENEGRAPH_INCLUDE_DIRS})
  else()
    message(FATAL_ERROR "OpenSceneGraph not found. Please set OSG_DIR.")
  endif()

  # Find osg
  find_package(osg)
  if(OSG_FOUND)
    set(MEPP_GUI_INC ${MEPP_GUI_INC} ${OSG_INCLUDE_DIR})
    set(MEPP_GUI_LIB ${MEPP_GUI_LIB} ${OSG_LIBRARIES})
  else()
    message(FATAL_ERROR "osg not found. Please set OSG_DIR.")
  endif()

  # Find osgViewer
  find_package(osgViewer)
  if(OSGVIEWER_FOUND)
    set(MEPP_GUI_INC ${MEPP_GUI_INC} ${OSGVIEWER_INCLUDE_DIR})
    set(MEPP_GUI_LIB ${MEPP_GUI_LIB} ${OSGVIEWER_LIBRARIES})
  else()
    message(FATAL_ERROR "osgViewer not found. Please set OSG_DIR.")
  endif()

  # Find osgUtil
  find_package(osgUtil)
  if(OSGUTIL_FOUND)
    set(MEPP_GUI_INC ${MEPP_GUI_INC} ${OSGUTIL_INCLUDE_DIR})
    set(MEPP_GUI_LIB ${MEPP_GUI_LIB} ${OSGUTIL_LIBRARIES})
  else()
    message(FATAL_ERROR "osgUtil not found. Please set OSG_DIR.")
  endif()

  # Find osgText
  find_package(osgText)
  if(OSGTEXT_FOUND)
    set(MEPP_GUI_INC ${MEPP_GUI_INC} ${OSGTEXT_INCLUDE_DIR})
    set(MEPP_GUI_LIB ${MEPP_GUI_LIB} ${OSGTEXT_LIBRARIES})
  else()
    message(FATAL_ERROR "osgText not found. Please set OSG_DIR.")
  endif()

  # Find osgGA
  find_package(osgGA)
  if(OSGGA_FOUND)
    set(MEPP_GUI_INC ${MEPP_GUI_INC} ${OSGGA_INCLUDE_DIR})
    set(MEPP_GUI_LIB ${MEPP_GUI_LIB} ${OSGGA_LIBRARIES})
  else()
    message(FATAL_ERROR "osgGA not found. Please set OSG_DIR.")
  endif()

  # Find osgDB
  find_package(osgDB)
  if(OSGDB_FOUND)
    set(MEPP_GUI_INC ${MEPP_GUI_INC} ${OSGDB_INCLUDE_DIR})
    set(MEPP_GUI_LIB ${MEPP_GUI_LIB} ${OSGDB_LIBRARIES})
  else()
    message(FATAL_ERROR "osgDB not found. Please set OSG_DIR.")
  endif()

  # Find osgFX
  find_package(osgFX)
  if(OSGFX_FOUND)
    set(MEPP_GUI_INC ${MEPP_GUI_INC} ${OSGFX_INCLUDE_DIR})
    set(MEPP_GUI_LIB ${MEPP_GUI_LIB} ${OSGFX_LIBRARIES})
  else()
    message(FATAL_ERROR "osgFX not found. Please set OSG_DIR.")
  endif()

  # Find osgShadow
  find_package(osgShadow)
  if(OSGSHADOW_FOUND)
    set(MEPP_GUI_INC ${MEPP_GUI_INC} ${OSGSHADOW_INCLUDE_DIR})
    set(MEPP_GUI_LIB ${MEPP_GUI_LIB} ${OSGSHADOW_LIBRARIES})
  else()
    message(FATAL_ERROR "osgShadow not found. Please set OSG_DIR.")
  endif()

  # Find osgWidget
  find_package(osgWidget)
  if(OSGWIDGET_FOUND)
    set(MEPP_GUI_INC ${MEPP_GUI_INC} ${OSGWIDGET_INCLUDE_DIR})
    set(MEPP_GUI_LIB ${MEPP_GUI_LIB} ${OSGWIDGET_LIBRARIES})
  else()
    message(FATAL_ERROR "osgWidget not found. Please set OSG_DIR.")
  endif()

  # Find osgManipulator
  find_package(osgManipulator)
  if(OSGMANIPULATOR_FOUND)
    set(MEPP_GUI_INC ${MEPP_GUI_INC} ${OSGMANIPULATOR_INCLUDE_DIR})
    set(MEPP_GUI_LIB ${MEPP_GUI_LIB} ${OSGMANIPULATOR_LIBRARIES})
  else()
    message(FATAL_ERROR "osgManipulator not found. Please set OSG_DIR.")
  endif()

  # Find osgQt - caution, no more 'osgQt' module embedded by default since openscenegraph 3.5.5 !
  message( STATUS "OSG version: " ${_osg_VERSION_MAJOR}.${_osg_VERSION_MINOR}.x )
  if( MSVC )
    set( OSGQT_INCLUDE_DIR
      "${PROJECT_SOURCE_DIR}/External/osgQt/osg34and36/include" )
    include_directories(${OSGQT_INCLUDE_DIR})
  elseif( APPLE )
    if( ${_osg_VERSION_MAJOR} EQUAL 3 AND ${_osg_VERSION_MINOR} EQUAL 5 )
      set( OSGQT_INCLUDE_DIR
        "${PROJECT_SOURCE_DIR}/External/osgQt/osg356/include" )
    else()
      set( OSGQT_INCLUDE_DIR
        "${PROJECT_SOURCE_DIR}/External/osgQt/osg34and36/include" )
    endif()
    include_directories(${OSGQT_INCLUDE_DIR})
  else() # Linux
    if( ${_osg_VERSION_MAJOR} EQUAL 3 AND ${_osg_VERSION_MINOR} LESS 4 )
      set( OSGQT_INCLUDE_DIR
        "${PROJECT_SOURCE_DIR}/External/osgQt/osg32/include" )
    elseif( ${_osg_VERSION_MAJOR} EQUAL 3 AND ${_osg_VERSION_MINOR} EQUAL 5 )
      set( OSGQT_INCLUDE_DIR
        "${PROJECT_SOURCE_DIR}/External/osgQt/osg356/include" )
    else()
      set( OSGQT_INCLUDE_DIR
        "${PROJECT_SOURCE_DIR}/External/osgQt/osg34and36/include" )
    endif()
    include_directories(${OSGQT_INCLUDE_DIR})
  endif()

  # Find OpenThreads # OpenThreads is part of OpenSceneGraph
  find_package(OpenThreads)
  if(OPENTHREADS_FOUND)
    set(MEPP_GUI_INC ${MEPP_GUI_INC} ${OPENTHREADS_INCLUDE_DIR})
    set(MEPP_GUI_LIB ${MEPP_GUI_LIB} ${OPENTHREADS_LIBRARY})
  else()
    message(FATAL_ERROR "OpenThreads not found. Please set OPENTHREADS_DIR.")
  endif()

  if( OSGVIEWER_FOUND )
    set(OPENSCENEGRAPH_FOUND_MEPP 1)
  else()
    message(FATAL_ERROR  "OpenSceneGraph not found.  Check the cmake variables associated to this package or disable it." )
  endif()

endif()

message(STATUS "--------------------------------------------------------------")
