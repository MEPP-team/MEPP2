if( NOT BUILD_USE_GUI OR NOT OPENSCENEGRAPH_FOUND_MEPP )
  return()
endif()

# Refer to DevelopersDesignCMakeFiles.dox entry 001
if( BUILD_USE_CGAL )
  include(${CGAL_USE_FILE})
else()
  include_directories( BEFORE "${PROJECT_SOURCE_DIR}/External/CGAL-4.14" )
endif()

if( BUILD_USE_OPENMESH )
  include_directories( ${OPENMESH_INCLUDE_DIR} )
endif()

if( BUILD_USE_PCL )
  include_directories( ${PCL_INCLUDE_DIRS} )
endif()

if( MSVC )
  # /FORCE:MULTIPLE for osgDB ifstream / ofstream conflict
  # -> here, really no choice... (because of OSG and OpenMesh already compiled libraries)
  set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /FORCE:MULTIPLE /IGNORE:4006 /IGNORE:4088" )
  #set( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /FORCE:MULTIPLE /IGNORE:4006 /IGNORE:4088" )
  # (Note that MSVC apparently considers 4088 too important to allow us to ignore it; I'm nevertheless leaving this here to document that we don't care about it.)

  # get rid of warning :
  # 'identifier' : decorated name length exceeded, name was truncated
  add_definitions( /wd4503 )

  # get rid (of most) of all random errors (occurs only with QT applications) :
  # 'Error 1 error MSB6006: "cmd.exe" exited with code 1'
  # -> /MP1 for ONLY ONE compilation unit (cl.exe); Without any value, it is set automatically depending on your number of main threads.
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP1")
endif()

set( UI_TESTS_SRC_DEMO
     mepp-gui )

set( GUILIB_DEMO ${MEPP_GUI_LIB} ${Boost_LIBRARIES} ${FBX_LIBRARY} )

### CGAL ###
if( BUILD_USE_CGAL )
  set( GUILIB_DEMO
    ${GUILIB_DEMO}
    ${CGAL_LIBRARY}
    ${CGAL_3RD_PARTY_LIBRARIES} )
endif()
### CGAL ###

### OPENMESH ###
if( BUILD_USE_OPENMESH )

  if(0)
    if( BUILD_USE_PCL ) ### For the moment, GUI needs a model of HalfedgeGraph. I (@jlevallois) have choose OpenMesh for the moment. This should be fixed.
      include_directories( ${PCL_INCLUDE_DIRS} )
      add_definitions( ${PCL_DEFINITIONS} )

      set( UI_TESTS_SRC_DEMO
         ${UI_TESTS_SRC_DEMO}
         testViewer_PointCloud )

      set( GUILIB_DEMO
        ${GUILIB_DEMO}
        ${PCL_COMMON_LIBRARIES} ${PCL_IO_LIBRARIES} )
    endif()
  endif()

    if( NOT BUILD_USE_CGAL )

      if( MSVC )
        # AUTOLINK is used for CGAL on Windows
        # Disable it for OpenMesh-only tests to prevent
        # false link error
        add_definitions( -DCGAL_NO_AUTOLINK_CGAL )
      endif()

      # We don't want our CGAL imports within External to constrain the non
      # CGAL implementations (e.g. OpenMesh) to be linked against CGAL libraries
      # just because of calls to CGAL's debugging/assertion... utilities (e.g.
      # CGAL::assertion_fail() ). Hence when in Debug mode the CGAL part of the
      # the Exertnals import is still instructed to be inhibited.
      # from being in debug mode.
      add_definitions( -DCGAL_NDEBUG )
    endif()

    set( GUILIB_DEMO
      ${GUILIB_DEMO}
      ${OPENMESH_LIBRARIES} )

endif()
### OPENMESH ###

### PCL ###
if( BUILD_USE_PCL )
  add_definitions( ${PCL_DEFINITIONS} )
  link_directories(${PCL_LIBRARY_DIRS})
  set( GUILIB_DEMO
       ${GUILIB_DEMO}
       ${PCL_LIBRARIES} )
endif()
### PCL ###

### VTK ###
if ( BUILD_USE_VTK )
  set( GUILIB_DEMO
       ${GUILIB_DEMO}
       ${VTK_LIBRARIES} )
endif()
### VTK ###

include_directories( "${MEPP_GUI_INC}" )
include_directories( "${CMAKE_CURRENT_BINARY_DIR}" )

# osgQt
if(MSVC)
  set( osgQt_SRC ${osgQt_SRC}
    "${PROJECT_SOURCE_DIR}/External/osgQt/osg34and36/GraphicsWindowQt.cxx" ) # see also /Cmake/GUIOptionalDependencies.cmake for include
elseif( APPLE )
  if( ${_osg_VERSION_MAJOR} EQUAL 3 AND ${_osg_VERSION_MINOR} EQUAL 5 )
    set( osgQt_SRC ${osgQt_SRC}
      "${PROJECT_SOURCE_DIR}/External/osgQt/osg356/GraphicsWindowQt.cxx" ) # see also /Cmake/GUIOptionalDependencies.cmake for include
  else()
    set( osgQt_SRC ${osgQt_SRC}
      "${PROJECT_SOURCE_DIR}/External/osgQt/osg34and36/GraphicsWindowQt.cxx" ) # see also /Cmake/GUIOptionalDependencies.cmake for include
  endif()
else() # Linux
  if( ${_osg_VERSION_MAJOR} EQUAL 3 AND ${_osg_VERSION_MINOR} LESS 4 )
    set( osgQt_SRC ${osgQt_SRC}
      "${PROJECT_SOURCE_DIR}/External/osgQt/osg32/GraphicsWindowQt.cxx" ) # see also /Cmake/GUIOptionalDependencies.cmake for include
  elseif( ${_osg_VERSION_MAJOR} EQUAL 3 AND ${_osg_VERSION_MINOR} EQUAL 5 )
    set( osgQt_SRC ${osgQt_SRC}
      "${PROJECT_SOURCE_DIR}/External/osgQt/osg356/GraphicsWindowQt.cxx" ) # see also /Cmake/GUIOptionalDependencies.cmake for include
  else()
    set( osgQt_SRC ${osgQt_SRC}
      "${PROJECT_SOURCE_DIR}/External/osgQt/osg34and36/GraphicsWindowQt.cxx" ) # see also /Cmake/GUIOptionalDependencies.cmake for include
  endif()
endif()

# common moc stuff to all executables... /
# IMPORTANT : automoc is off for qt5 because we wan't same behavior
# with qt4 and/or qt5...
set( moc_files ${moc_files}
  "${PROJECT_SOURCE_DIR}/Visualization/Plugins/PluginDialog.hpp"
#  "${PROJECT_SOURCE_DIR}/Visualization/Qt/SimpleLineEdit.h"
#  "${PROJECT_SOURCE_DIR}/Visualization/Qt/SimpleCheckBox.h"
  "${PROJECT_SOURCE_DIR}/Visualization/SimpleWindow.h"
  )
#set( moc_files ${moc_files} "add_here_your_header_with_Q_OBJECT.h" )

# ui
set( ui_files ${ui_files}
  "${PROJECT_SOURCE_DIR}/Visualization/SimpleWindow.ui" )
#set( ui_files ${ui_files} "add_here_your_ui_file" )

# qrc
set( qrc_files ${qrc_files}
  "${PROJECT_SOURCE_DIR}/Visualization/Applications/mepp-gui.qrc")
#set( qrc_files ${qrc_files} "add_here_your_qrc_file" )

if( BUILD_USE_QT5 )
  qt5_wrap_cpp( MOC_FILES ${moc_files} )
  qt5_wrap_ui( UI_FILES ${ui_files} )
  qt5_add_resources( QT_RESOURCE_FILES ${qrc_files} )
else()
  qt4_wrap_cpp( MOC_FILES ${moc_files}
    OPTIONS -DBOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
    -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED )
  qt4_wrap_ui( UI_FILES ${ui_files} )
  qt4_add_resources( QT_RESOURCE_FILES ${qrc_files} )
endif()

# For Windows set the .ico file with the .rc file
if( MSVC )
  set( resource_files ${resource_files}
    "${PROJECT_SOURCE_DIR}/Visualization/Applications/resources/mepp-gui.rc" )
endif()

# For Apple set the .icns file containing icons
if(APPLE)
  # set how it shows up in the Info.plist file
  set(MACOSX_BUNDLE_ICON_FILE mepp-gui.icns)

  # set where in the bundle to put the icns file
  set_source_files_properties(
    "${PROJECT_SOURCE_DIR}/Visualization/Applications/resources/mepp-gui.icns"
    PROPERTIES MACOSX_PACKAGE_LOCATION Resources)

  # include the icns file in the target
  set( resource_files ${resource_files}
    "${PROJECT_SOURCE_DIR}/Visualization/Applications/resources/mepp-gui.icns" )
endif()

# debug visualization in CI
# uncomment to enable
#if( (DEFINED ENV{GITLAB_CI}) OR (DEFINED ENV{TRAVIS}) OR (DEFINED ENV{APPVEYOR}) )
#  add_definitions( -DDEBUG_VISU2 )
#endif()

foreach(FILE ${UI_TESTS_SRC_DEMO})

  add_executable( ${FILE} MACOSX_BUNDLE
    ${osgQt_SRC}
    ${MOC_FILES}
    ${QT_RESOURCE_FILES}
    ${resource_files}
    ${FILE}.cpp )

  if("${FILE}" STREQUAL "mepp-gui")
    target_sources(mepp-gui PUBLIC ${UI_FILES})
  else()
    add_dependencies(${FILE} mepp-gui)
  endif()

  target_link_libraries( ${FILE} ${GUILIB_DEMO} )

  if("${FILE}" STREQUAL "testViewer_PointCloud")
    set( DATA_TEST
        ${TESTING_DATA_DIR}/ism_test_cat.pcd false)
  else() # not PCL
    set( DATA_TEST
        ${TESTING_DATA_DIR}/casting.off )
  endif()

  if( APPLE )
    if( (NOT DEFINED ENV{TRAVIS}) AND (NOT DEFINED ENV{GITLAB_CI}) ) # OS X - random freeze under travis and GITLAB_CI (but only under OSX VMs)
      add_test( ${FILE}
                ${FILE}.app/Contents/MacOS/${FILE}
                ${DATA_TEST}
                test )
    endif()
  elseif( MSVC )
    add_test( ${FILE}
              ${FILE}
              ${DATA_TEST}
              ${VTEST} )
  else()
    # LINUX
    add_test( ${FILE}
              ${FILE}
              ${DATA_TEST}
              test )
    # display backtrace in case of SEGFAULT ; GDB required
    # uncomment to activate
    #add_test( ${FILE}
    #          "${PROJECT_SOURCE_DIR}/Tools/segfault_debugger.sh"
    #          "./${FILE}"
    #          ${DATA_TEST}
    #          test )
  endif()

endforeach()

# --- EXECUTABLE(S) ---

if( MSVC )
  if( QT4_FOUND_MEPP )
    set( QT_PLUGINS_DIR "${QTDIR}/plugins" )
  elseif( QT5_FOUND_MEPP )
    set( QT_PLUGINS_DIR "${QT5_DIR}/plugins" )
  endif()

  #---------------------------------------------------------------------------
  # Copy Qt plugins to 'Debug & Release' directories & configure qt.conf file
  #---------------------------------------------------------------------------
  file( GLOB qtplugin_dirs RELATIVE
    "${QT_PLUGINS_DIR}"
    "${QT_PLUGINS_DIR}/imageformats*"
    "${QT_PLUGINS_DIR}/platforms*" )
  file( REMOVE_RECURSE
    "${CMAKE_CURRENT_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/Debug/QtPlugins" )
  file( REMOVE_RECURSE
    "${CMAKE_CURRENT_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/Release/QtPlugins" )
  file( MAKE_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/Debug/QtPlugins" )
  file( MAKE_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/Release/QtPlugins" )
  foreach( qtplugin ${qtplugin_dirs} )
    file( COPY "${QT_PLUGINS_DIR}/${qtplugin}"
      DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/Debug/QtPlugins" )
    file( COPY "${QT_PLUGINS_DIR}/${qtplugin}"
      DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/Release/QtPlugins" )
  endforeach()
  file( WRITE
    "${CMAKE_CURRENT_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/Debug/qt.conf"
    "[Paths]\nPlugins = QtPlugins" )
  file( WRITE
    "${CMAKE_CURRENT_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/Release/qt.conf"
    "[Paths]\nPlugins = QtPlugins" )
endif()

#---------------------------------------------------------------------------
# Copying shaders
#---------------------------------------------------------------------------
set( SHADERS_DIR "${PROJECT_SOURCE_DIR}/Visualization/Shaders" )
file( GLOB SHADERS_FILES RELATIVE
  "${SHADERS_DIR}"
  "${SHADERS_DIR}/*" )
if( MSVC )
	file( MAKE_DIRECTORY
	  "${CMAKE_CURRENT_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/Debug/Shaders" )
	foreach( shader ${SHADERS_FILES} )
	  file( COPY "${SHADERS_DIR}/${shader}"
	    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/Debug/Shaders" )
	endforeach()
	file( MAKE_DIRECTORY
	  "${CMAKE_CURRENT_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/Release/Shaders" )
	foreach( shader ${SHADERS_FILES} )
	  file( COPY "${SHADERS_DIR}/${shader}"
	    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/Release/Shaders" )
	endforeach()
else()
	file( MAKE_DIRECTORY
	  "${CMAKE_CURRENT_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/Shaders" )
	foreach( shader ${SHADERS_FILES} )
	  file( COPY "${SHADERS_DIR}/${shader}"
	    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/Shaders" )
	endforeach()
endif()

# --- Visualization PLUGINS(S) ---

# Set location of plugins : IMPORTANT !
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})
# Set location of plugins : IMPORTANT !

include( "${CMAKE_CURRENT_SOURCE_DIR}/../PluginFilters/PluginFilters.cmake" )

#####
# display CMake include directories for debugging purpose
#get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
#foreach(dir ${dirs})
#  message(STATUS "dir='${dir}'")
#endforeach()


##TODO-elo make a function
## variables to be exported to external projects

# include dirs
get_property( INCLUDE_DIRECTORIES_TMP
              DIRECTORY
              PROPERTY INCLUDE_DIRECTORIES
              )
set( MEPP2_GUI_INCLUDE_DIRECTORIES
     "${INCLUDE_DIRECTORIES_TMP}"
     CACHE INTERNAL "no-desc" FORCE
     )

# system include dirs
get_property( SYSTEM_INCLUDE_DIRECTORIES_TMP
              DIRECTORY
              PROPERTY SYSTEM_INCLUDE_DIRECTORIES
              )
set( MEPP2_GUI_SYSTEM_INCLUDE_DIRECTORIES
     "${SYSTEM_INCLUDE_DIRECTORIES_TMP}"
     CACHE INTERNAL "no-desc" FORCE
     )

# libraries to link with
set( MEPP2_GUI_LIBRARIES
     ${GUILIB_DEMO}
     CACHE INTERNAL "no-desc" FORCE
     )

# link libraries (to link all target with, silently!)
get_property( LINK_LIBRARIES_TMP
              DIRECTORY
              PROPERTY LINK_LIBRARIES
              )
set( MEPP2_GUI_LINK_LIBRARIES
     "${LINK_LIBRARIES_TMP}"
     CACHE INTERNAL "no-desc" FORCE
     )

# link dirs
get_property( LINK_DIRECTORIES_TMP
              DIRECTORY
              PROPERTY LINK_DIRECTORIES
              )
set( MEPP2_GUI_LINK_DIRECTORIES
     "${LINK_DIRECTORIES_TMP}"
     CACHE INTERNAL "no-desc" FORCE
     )

# compilation flags and define
get_property( COMPILE_DEFINITIONS_TMP
              DIRECTORY
              PROPERTY COMPILE_DEFINITIONS
              )
set( MEPP2_GUI_COMPILE_DEFINITIONS
     ${COMPILE_DEFINITIONS_TMP}
     CACHE INTERNAL "no-desc" FORCE
     )

## extra variables to be exported to external projects
## specific to plugins

set( MEPP2_osgQt_SRC
     ${osgQt_SRC}
     CACHE INTERNAL "no-desc" FORCE
     )
