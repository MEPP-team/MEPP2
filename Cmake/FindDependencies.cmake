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
  if( ${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION} VERSION_LESS "4.14" )
    message (FATAL_ERROR "CGAL required minimum version is 4.14 - Found: ${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION}" )
    return()
  endif()
  add_definitions( -DFEVV_USE_CGAL )
  # FIXME: used in if and else (see below). Factorize this !
  add_definitions( -DBOOST_ALL_DYN_LINK )
else()
  add_definitions( -DBOOST_ALL_DYN_LINK )
  # Refer to Doc/Devel/CMakeFiles.md entry 003
  add_definitions( -DCGAL_NO_AUTOLINK_CGAL )
  add_definitions( -DCGAL_NDEBUG )
  add_definitions( -DCGAL_HEADER_ONLY )
  add_definitions( -DCGAL_DISABLE_GMP )
  message ( STATUS "Use local LGPL CGAL")
endif()

##### Boost package finding (mandatory)
find_package(Boost COMPONENTS thread system filesystem REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
else()
  message(FATAL_ERROR "Boost not found. Please install it.")
endif(Boost_FOUND)

##### EIGEN3 package finding (mandatory)
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

##### IMG-3RDPARTY libraries finding
if( BUILD_USE_IMG-3RDPARTY )
  FIND_PACKAGE(JPEG)
  if ( JPEG_FOUND AND EXISTS ${JPEG_LIBRARIES} )
    set(FEVV_HAS_ONE_IMG_LIBRARY 1)

    add_definitions( -DFEVV_USE_JPEG )
    include_directories( ${JPEG_INCLUDE_DIR} )

    set( IMG-3RDPARTY_LIBRARIES ${IMG-3RDPARTY_LIBRARIES} ${JPEG_LIBRARIES} )

    #message( STATUS "---> JPEG_INCLUDE_DIR: ${JPEG_INCLUDE_DIR}" )
    #message( STATUS "---> JPEG_LIBRARIES: ${JPEG_LIBRARIES}" )
  else()
    #message ( "Unfound Jpeg library." )
  endif ()

  FIND_PACKAGE(PNG) # and ZLIB
  if ( PNG_FOUND )
	  if ( ZLIB_FOUND AND EXISTS ${ZLIB_LIBRARIES} )
	    set(FEVV_HAS_ONE_IMG_LIBRARY 1)

	    add_definitions( -DFEVV_USE_PNG )
	    include_directories( ${PNG_PNG_INCLUDE_DIR} )

	    set( IMG-3RDPARTY_LIBRARIES ${IMG-3RDPARTY_LIBRARIES} ${PNG_LIBRARIES} )

	    #message( STATUS "---> PNG_PNG_INCLUDE_DIR: ${PNG_PNG_INCLUDE_DIR}" )
	    #message( STATUS "---> PNG_LIBRARIES: ${PNG_LIBRARIES}" )
	  else()
	    #message ( "Unfound Png library." )
	  endif ()
  endif ()

  FIND_PACKAGE(TIFF)
  if ( TIFF_FOUND AND EXISTS "${TIFF_LIBRARIES}" )
    set(FEVV_HAS_ONE_IMG_LIBRARY 1)

    add_definitions( -DFEVV_USE_TIFF )
    include_directories( ${TIFF_INCLUDE_DIR} )

    set( IMG-3RDPARTY_LIBRARIES ${IMG-3RDPARTY_LIBRARIES} ${TIFF_LIBRARIES} )

    #message( STATUS "---> TIFF_INCLUDE_DIR: ${TIFF_INCLUDE_DIR}" )
    #message( STATUS "---> TIFF_LIBRARIES: ${TIFF_LIBRARIES}" )
  else()
    #message ( "Unfound Tiff library." )
  endif ()

  if ( NOT FEVV_HAS_ONE_IMG_LIBRARY )
    message (FATAL_ERROR "None of Jpeg or Png (or Zlib) or Tiff library found. Turn BUILD_USE_IMG-3RDPARTY to OFF.")
  else ()
    if ( NOT MSVC )
        add_definitions( -Dcimg_display=0 ) # Flags used to disable display capablities of CImg (for Linux and Mac OS)
    endif ()
    if ( DEFINED ENV{TRAVIS} OR DEFINED ENV{APPVEYOR} )
        add_definitions( -Dcimg_display=0 ) # Flags used to disable display capablities of CImg (for all CIs)
    endif ()
  endif ()
endif()

##### PCL package finding
if( BUILD_USE_PCL )
  find_package(PCL 1.3 REQUIRED COMPONENTS common io search features)
  if ( PCL_FOUND )
    add_definitions( -DFEVV_USE_PCL )
  else()
    message (FATAL_ERROR "PCL not found. Turn BUILD_USE_PCL to OFF.")
  endif ()

  # not necessary for PCL usage, just here for flann example...
  find_package(Flann)
  if ( FLANN_FOUND )
    add_definitions( -DFEVV_USE_FLANN )
  else()
    message ( "Unfound Flann package." )
  endif ()
endif ()

##### VTK package finding
if( BUILD_USE_VTK )
  find_package(VTK REQUIRED)
  if ( VTK_FOUND )
    include_directories( ${VTK_INCLUDE_DIRS} )
    add_definitions( -DFEVV_USE_VTK )
  else()
    message (FATAL_ERROR "VTK not found. Turn BUILD_USE_VTK to OFF.")
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
      set( FBX_LIBRARY ${FBX_LIBRARY}/${BUILD_CONFIG}/libfbxsdk.so ${CMAKE_DL_LIBS} )
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

##### Draco package finding
if( BUILD_USE_DRACO )
  find_package(Draco)
  if ( draco_FOUND )
    include_directories( ${draco_INCLUDE_DIRS} )
    add_definitions( -DFEVV_USE_DRACO )
  else()
    message( FATAL_ERROR
             "Draco not found. Please set DRACO_DIR or turn BUILD_USE_DRACO to OFF.")
  endif ()
endif ()
