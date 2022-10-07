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

# -----------------------------------------------------------------------------
# BUILD_USE_QT6:
#
# 1)  For the moment, under Ubuntu 22.04, VTK must be set to OFF if BUILD_USE_QT6 is ON (because VTK is linked with Qt5 !)
# 2a) For the moment, under Windows*, the BUILD_USE_QT6 option is hidden because of additional dependencies for Eigen and OpenMesh
#     (newer versions: Eigen 3.3.9 and OpenMesh 8.1)
# 2b) For the moment, under Windows* (current binary kit), PCL must be set to OFF if BUILD_USE_QT6 is ON (because Flann 1.9.1 is too old !)
#
#     * BUILD_USE_QT6 under Windows requires MSVC 2019 !
# -----------------------------------------------------------------------------
#option( BUILD_USE_QT6 "Using Qt6." OFF ) # TODO, comment

if( BUILD_USE_QT6 )
  set( BUILD_USE_QT5 OFF )
  set( LIST_OPTION ${LIST_OPTION} [QT6]\ ) # ??
  message(STATUS "   BUILD_USE_QT6  true   (Using of Qt6 instead of Qt4)")
else()
  message(STATUS "   BUILD_USE_QT6  false  (Using of Qt6 instead of Qt4)")
endif()

if( BUILD_USE_QT5 )
  set( BUILD_USE_QT6 OFF )
  set( LIST_OPTION ${LIST_OPTION} [QT5]\ ) # ??
  message(STATUS "   BUILD_USE_QT5  true   (Using of Qt5 instead of Qt4)")
else()
  message(STATUS "   BUILD_USE_QT5  false  (Using of Qt5 instead of Qt4)")
endif()

if( BUILD_USE_GUI )
  set( LIST_OPTION ${LIST_OPTION} [UI]\ ) # ??
  message(STATUS "   BUILD_USE_GUI  true   (Build user interface)")
else()
  message(STATUS "   BUILD_USE_GUI  false  (Don't build user interface)")
endif()

# -----------------------------------------------------------------------------
# Look for Qt (needed by openscenegraph visualization).
# -----------------------------------------------------------------------------
set( QT4_FOUND_MEPP 0 )
set( QT5_FOUND_MEPP 0 )
set( QT6_FOUND_MEPP 0 )
if( BUILD_USE_GUI )
  if( BUILD_USE_QT6 )
    # QT6 Handling
    if(DEFINED ENV{QT6_DIR})
      set( QT6_DIR $ENV{QT6_DIR} )
    endif()
    set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${QT6_DIR})

    # Instruct CMake to NOT run moc automatically when needed
    set(CMAKE_AUTOMOC OFF)

    find_package(Qt6 COMPONENTS Widgets OpenGL OpenGLWidgets Xml REQUIRED) # for Qt6 add 'OpenGLWidgets'

    if( Qt6Widgets_FOUND AND Qt6OpenGL_FOUND AND Qt6OpenGLWidgets_FOUND AND Qt6Xml_FOUND ) # for Qt6 add 'AND Qt6OpenGLWidgets_FOUND'
      set(QT6_FOUND_MEPP 1)
      message( STATUS "Qt6 (Widgets, OpenGL, OpenGLWidgets, Xml modules) found (needed by OpenSceneGraph compiled with Qt6)." ) # for Qt6 add ', OpenGLWidgets'

      add_definitions("-DFEVV_USE_QT5") # important !!!
      add_definitions("-DFEVV_USE_QT6") # for future... ???

      set(MEPP_GUI_LIB ${MEPP_GUI_LIB}
        ${Qt6Widgets_LIBRARIES}
        ${Qt6OpenGL_LIBRARIES}
        ${Qt6OpenGLWidgets_LIBRARIES} # for Qt6
        ${Qt6Xml_LIBRARIES})
      set(MEPP_GUI_INC ${MEPP_GUI_INC}
        ${Qt6Widgets_INCLUDES_DIRS}
        ${Qt6OpenGL_INCLUDES_DIR}
        ${Qt6OpenGLWidgets_INCLUDES_DIR} # for Qt6
        ${Qt6Xml_INCLUDES_DIR})

    else()
      message( STATUS "One of Qt6's modules was not found (needed by OpenSceneGraph)." )
    endif()
  elseif( BUILD_USE_QT5 ) # sample for help : QT5_DIR=/usr/local/Cellar/qt/5.9.2 or Qt5_DIR=/usr/local/Cellar/qt/5.9.2/lib/cmake/Qt5 (caution, here 't' not 'T' !)
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
      message( STATUS "Qt5 (Widgets, OpenGL, Xml modules) found (needed by OpenSceneGraph compiled with Qt5)." )

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
  # --> With Qt5 and Qt6, now, the new module is 'osgQOpenGL'
  message("--> OSG version: " ${_osg_VERSION_MAJOR}.${_osg_VERSION_MINOR}.x )
  if( MSVC )
    if( BUILD_USE_QT5 OR BUILD_USE_QT6 ) # FOR_QT6
      set( OSGQT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/External/osgQOpenGL/osg34and36/include" )
    else()
      set( OSGQT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/External/osgQt/osg34and36/include" )
    endif()
    include_directories(${OSGQT_INCLUDE_DIR})
  elseif( APPLE )
    if( BUILD_USE_QT5 OR BUILD_USE_QT6 ) # FOR_QT6
      set( OSGQT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/External/osgQOpenGL/osg34and36/include" )
    else()
      if( ${_osg_VERSION_MAJOR} EQUAL 3 AND ${_osg_VERSION_MINOR} EQUAL 5 )
        set( OSGQT_INCLUDE_DIR
          "${PROJECT_SOURCE_DIR}/External/osgQt/osg356/include" )
      else()
        set( OSGQT_INCLUDE_DIR
          "${PROJECT_SOURCE_DIR}/External/osgQt/osg34and36/include" )
      endif()
    endif()
    include_directories(${OSGQT_INCLUDE_DIR})
  else() # Linux
    if( BUILD_USE_QT5 OR BUILD_USE_QT6 ) # FOR_QT6
      set( OSGQT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/External/osgQOpenGL/osg34and36/include" )
    else()
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
