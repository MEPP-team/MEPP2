##### CGAL package finding
# Caveat emptor: for some undocumented reason CGAL package detection messes
#     up (it clears it completely at least on OSX) the content of the
#     Boost_LIBRARIES variable as set up by the Boost package detection.
#     In order to avoid this side effect, CGAL package detection must be
#     placed BEFORE Boost package detection.
if( BUILD_USE_CGAL )
  add_definitions( -DCGAL_EIGEN3_ENABLED ) ##### SOME CGAL FUNCTIONALITIES DEPEND ON EIGEN3

  set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "Prevent CGAL to override CMAKE flags")
  # The find_package( CGAL 4.7 REQUIRED ) directive fails on:
  #  * OSX 10.11 with a homebrew install
  #  * Ubuntu (as offered by Travis).
  # This might be due to CGAL config file that doesn't seem to positionate
  # version numbers properly. Indeed cmake 3.4 expects the config file to define
  # variable like "<package>_VERSION_MAJOR" whereas (at least on OSX 10.11 with
  # Homebrew packages) CGALConfig.cmake defines CGAL_MAJOR_VERSION
  # Refer to https://cmake.org/cmake/help/v3.4/command/find_package.html
  # Hence the following manual check
  find_package( CGAL REQUIRED CONFIG )
  if( ${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION} VERSION_LESS "4.11" )
    message ( "CGAL required minimum version is 4.11 - Found: ${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION}" )
    return()
  endif()
  add_definitions( -DFEVV_USE_CGAL )
  # FIXME: used in if and else (see below). Factorize this !
  add_definitions( -DBOOST_ALL_DYN_LINK )

  # CGAL 4.11.0 requires to be patched (4.11.1 should fix this trouble)
  if( NOT WIN32 )
    if( ${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION}.${CGAL_BUGFIX_VERSION}
        VERSION_EQUAL "4.11.0" )
      file(MD5 ${CGAL_INCLUDE_DIRS}/CGAL/boost/graph/io.h cgal_io_h_md5)
      if( "${cgal_io_h_md5}" STREQUAL "58b04f1e2334d3136276469531372eed")
        message(FATAL_ERROR "CGAL 4.11.0 must be patched: refer to Install.md.")
      endif()
    endif()
  endif()
else()
  add_definitions( -DBOOST_ALL_DYN_LINK )
  # Refer to Doc/Devel/CMakeFiles.md entry 003
  add_definitions( -DCGAL_NO_AUTOLINK_CGAL )
  add_definitions( -DCGAL_NDEBUG )
  message ( STATUS "Use local LGPL CGAL")
endif()

##### Boost is a must unless when building the documentation only
if( BUILD_USE_CGAL OR BUILD_USE_OPENMESH OR BUILD_USE_AIF )
  find_package(Boost COMPONENTS thread system filesystem REQUIRED)
  if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
  endif(Boost_FOUND)
endif()

##### EIGEN3
if (WIN32)
	if(DEFINED ENV{EIGEN3_INCLUDE_DIR})
		SET( EIGEN3_INCLUDE_DIR $ENV{EIGEN3_INCLUDE_DIR} )
	ENDIF()
endif (WIN32)
# Eigen is needed except for the documentation
if( BUILD_USE_CGAL OR BUILD_USE_OPENMESH OR BUILD_USE_AIF )
	FIND_PACKAGE(Eigen3)
	if(EIGEN3_FOUND)
	  INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
	ELSE(EIGEN3_FOUND)
	  MESSAGE(FATAL_ERROR
			  "Eigen3 not found. Please set EIGEN3_INCLUDE_DIR to the root directory of your Eigen3 installation.")
	endif(EIGEN3_FOUND)
endif()

##### OpenMesh package finding
if( BUILD_USE_OPENMESH )
  find_package( OpenMesh REQUIRED )
  if ( OPENMESH_FOUND )
    add_definitions( -DFEVV_USE_OPENMESH )
  else()
    message( FATAL_ERROR
             "OpenMesh not found. Please set OPENMESH_DIR or turn BUILD_USE_OPENMESH to OFF.")
  endif ()
endif ()

##### AIF
if( BUILD_USE_AIF )
  add_definitions( -DFEVV_USE_AIF )
endif()

##### PCL package finding
if( BUILD_USE_PCL )
  find_package(PCL 1.3 REQUIRED COMPONENTS common io)
  if ( PCL_FOUND )
    add_definitions( -DFEVV_USE_PCL )
  else()
    message ( "Unfound PCL package: try with PCL support disabled.")
    message ( "(try toggling BUILD_USE_PCL to OFF)." )
  endif ()

  # not necessary for PCL usage, just here for flann example...
  find_package(Flann)
  if ( FLANN_FOUND )
    add_definitions( -DFEVV_USE_FLANN )
  else()
    message ( "Unfound Flann package.")
  endif ()
endif ()

##### VTK package finding
if( BUILD_USE_VTK )
  find_package(VTK REQUIRED)
  if ( VTK_FOUND )
    include_directories( ${VTK_INCLUDE_DIRS} )
    add_definitions( -DFEVV_USE_VTK )
  else()
    message ( "Unfound VTK package: try with VTK support disabled.")
    message ( "(try toggling BUILD_USE_VTK to OFF)." )
  endif ()
endif ()

##### FBX package finding - primary support, not finished ! (use /Cmake/FindFBX.cmake instead after having done 3x 'Os updates & tests' with it...)
if( BUILD_USE_FBX )
  if( DEFINED ENV{FBX_DIR} )
    set( FBX_DIR $ENV{FBX_DIR} )
  endif()
  if( EXISTS ${FBX_DIR} )
    set( FBX_INCLUDE_DIRS ${FBX_DIR}/include )
    include_directories( ${FBX_INCLUDE_DIRS} )

    set( FBX_LIBRARY       ${FBX_DIR}/lib )

    if( "${CMAKE_BUILD_TYPE}" MATCHES "(D|d)eb" )
      set ( BUILD_CONFIG "debug" )
    elseif( "${CMAKE_BUILD_TYPE}" MATCHES "(R|r)el" )
      set ( BUILD_CONFIG "release" )
    else()
      message( "Unrecognized build type to use FBX." )
    endif()

    if(UNIX AND NOT APPLE) # Linux
      set( FBX_LIBRARY ${FBX_LIBRARY}/${BUILD_CONFIG}/libfbxsdk.so )
    elseif(APPLE) # Mac OS
      set( FBX_LIBRARY ${FBX_LIBRARY}/${BUILD_CONFIG}/libfbxsdk.dylib )
    else() # Windows
      set( FBX_LIBRARY ${FBX_LIBRARY}/${BUILD_CONFIG}/libfbxsdk.lib )
    endif()

    add_definitions( -DFEVV_USE_FBX )
    add_definitions( -DFBXSDK_SHARED )

    message( STATUS "FBX SDK loaded." )
  else()
    message( "FBX SDK not found." )
  endif()
endif()
